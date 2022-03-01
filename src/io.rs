use std::io::Error as IoError;
use std::path::Path;

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use noodles_bam::record::Record as NoodlesRecord;
use noodles_bam::{AsyncReader as NoodlesAsyncReader, AsyncWriter as NoodlesAsyncWriter};
use noodles_bgzf::AsyncReader as BgzfAsyncReader;
use noodles_bgzf::AsyncWriter as BgzfAsyncWriter;
use noodles_sam::header::{self, Header as BamHeader, ReferenceSequences};
use thiserror::Error;
use tokio::fs::File;
use tokio::io::{self, AsyncRead, AsyncWrite};

use crate::record::{FragmentCoord, UmiRecord, Position};
use crate::bktree::{BkTree, Dist};
use crate::optical::*;

pub struct BamIo<R, W>
where
    R: AsyncRead,
    W: AsyncWrite,
{
    in_bam: NoodlesReader<R>,
    out_bam: NoodlesWriter<W>,
    next_record: Option<NoodlesRecord>,
    reference_sequences: ReferenceSequences,
}

type NoodlesReader<R> = NoodlesAsyncReader<BgzfAsyncReader<R>>;
type NoodlesWriter<W> = NoodlesAsyncWriter<BgzfAsyncWriter<W>>;

/// BamIo is created from a `Read` and A `Write`. Upon construction it reads a BAM header and
/// reference sequence from the reader and writes the modified information to the writer. BamIo can
/// then reads records or bundles or records from the provided BAM file reader. The bam file must
/// be sorted to correctly read bundles of reads at the same position
impl<R, W> BamIo<R, W>
where
    R: AsyncRead + std::marker::Unpin,
    W: AsyncWrite + std::marker::Unpin,
{
    /// Create a new reader for `reader`. `reader` should be a at the beginning of a BAM file. The
    /// BAM header will be read during construction.
    pub async fn new(read: R, write: W) -> Result<BamIo<R, W>, BamIoError> {

        let mut in_bam = NoodlesAsyncReader::new(read);
        let header = in_bam.read_header().await?;
        let reference_sequences = in_bam.read_reference_sequences().await?;

        let mut out_bam = out_bam(write, &header, &reference_sequences).await?;

        Ok(BamIo {
            in_bam,
            out_bam,
            next_record: None,
            reference_sequences,
        })
    }


    /// Return the next read from the bam read stream or the peeked value
    pub async fn read(&mut self) -> io::Result<Option<NoodlesRecord>> {
        if let Some(record) = self.next_record.take() {
            Ok(Some(record))
        } else {
            self.read_from_bam().await
        }
    }

    /// Return the next read from the bam read stream or the peeked value if
    /// the position matches the provided position
    pub async fn read_at_pos(&mut self, pos: Position) -> io::Result<Option<NoodlesRecord>> {
        if self.next_record.is_none() {
            self.next_record = self.read_from_bam().await?;
        }

        let next_pos = self.next_record.as_ref().unwrap().position();
        if next_pos == Some(pos) {
            Ok(Some(self.next_record.take().unwrap()))
        } else {
            Ok(None)
        }
    }

    /// Read a new record from the stream and return it. if A record was previously read into
    /// the peek location this record is discarded and re-used.
    async fn read_from_bam(&mut self) -> io::Result<Option<NoodlesRecord>> {
        let mut record = NoodlesRecord::default();
        match self.in_bam.read_record(&mut record).await? {
            0 => Ok(None),
            _n => Ok(Some(record)),
        }
    }

    pub async fn read_bundle(&mut self, bundle: &mut Vec<NoodlesRecord>) -> io::Result<bool> {
        bundle.clear();
        if let Some(first) = self.read().await? {
            let pos = first.position();
            bundle.push(first);
            if let Some(pos) = pos {
                while let Some(next) = self.read_at_pos(pos).await? {
                    bundle.push(next);
                }
            }
        }
        Ok(!bundle.is_empty())
    }

    pub async fn read_bundle_owned(&mut self) -> io::Result<Vec<NoodlesRecord>> {
        let mut bundle = Vec::new();
        if let Some(first) = self.read().await? {
            let pos = first.position();
            bundle.push(first);
            if let Some(pos) = pos {
                while let Some(next) = self.read_at_pos(pos).await? {
                    bundle.push(next);
                }
            }
        }
        Ok(bundle)
    }

    pub async fn write_record(&mut self, record: &NoodlesRecord) -> io::Result<()> {
        self.out_bam.write_record(record).await
    }

    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.out_bam.shutdown().await
    }
}

async fn out_bam<W: AsyncWrite + std::marker::Unpin>(w: W, header: &str, reference_sequences: &ReferenceSequences, ) -> Result<NoodlesWriter<W>, BamIoError> {
    let mut writer = NoodlesAsyncWriter::new(w);
    let mut header: BamHeader = header.parse()
       .map_err(|e: header::ParseError| BamIoError::ParseError(e.to_string()))?;
    header
        .programs_mut()
        .insert("umidedup".to_string(), header::Program::new("umidedup"));
    writer.write_header(&header).await?;
    writer
        .write_reference_sequences(reference_sequences)
        .await?;

    Ok(writer)
}

#[derive(Debug, Error)]
pub enum BamIoError {
    #[error("Io error reading BAM")]
    IoError(#[from] std::io::Error),
    #[error("ParseError reading BAM: {0}")]
    ParseError(String),
    //#[error("Error using record")]
    //RecordError(#[from] RecordError),
}

#[derive(Debug, Eq, PartialEq, Hash)]
pub enum FragmentStart {
    Fw(i32),
    Rev(i32),
    TemplateEnd(i32),
}

#[derive(Debug, Eq, PartialEq, Hash)]
pub enum FragmentEnd {
    Fw(i32),
    Rev(i32),
    TemplateEnd(i32),
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
struct UmiIndex {
    index: usize,
    umi: Vec<u8>,
}

impl Dist for UmiIndex {
    fn dist(&self, b: &UmiIndex) -> usize {
        let la = self.umi.len();
        let lb = b.umi.len();
        if la == lb {
            self.umi
                .iter()
                .zip(b.umi.iter())
                .filter(|&(a, b)| a == &b'N' || b == &b'N' || a != b)
                .count()
        } else {
            std::cmp::max(la, lb)
        }
    }
}

#[derive(Debug, Default)]
pub struct ReadPartitions {
    pub seen: AHashMap<Vec<u8>, MarkResult>,
    pub current_reads: Vec<UmiRecord>,
    partitions: AHashMap<(FragmentCoord, FragmentCoord), BkTree<UmiIndex>>,
    // Mates wait for their mate to be decided.
    in_partitions: AHashSet<Vec<u8>>,
    current_read_mates: Vec<UmiRecord>,
    pub optical_dup_count: usize,
}

#[derive(Debug)]
pub enum RecordStatus {
    Unusable(UmiRecord),
    PossibleDup,
    Duplicate,
    SeenMate(UmiRecord),
    InPartition,
}

#[derive(Debug)]
pub struct MarkResult {
    corrected_umi: Option<Vec<u8>>,
    is_dup: bool,
}

impl ReadPartitions {
    pub fn add_record(&mut self, mut r: UmiRecord) -> RecordStatus {
        if r.is_paired() {
            if let Some(result) = self.seen.remove(r.read_name()) {
                if result.is_dup {
                    r.flag_dup();
                }
                if let Some(umi) = result.corrected_umi {
                    r.correct_barcode_tag(String::from_utf8_lossy(&umi));
                }
                return RecordStatus::SeenMate(r);
            } else if self.in_partitions.contains(r.read_name()) {
                self.current_read_mates.push(r);
                return RecordStatus::InPartition;
            }
        }

        //determine fragment coordinates
        if let Some(fragment) = r.fragment_markers() {
            let record_idx = self.current_reads.len();
            let umi = UmiIndex {
                umi: r.umi.to_vec(),
                index: record_idx,
            };
            self.in_partitions.insert(r.read_name().to_owned());
            self.current_reads.push(r);

            //insert into tree
            let e = self.partitions.entry(fragment).or_insert_with(BkTree::new);
            if e.insert(umi) {
                RecordStatus::PossibleDup
            } else {
                RecordStatus::Duplicate
            }
        } else {
            RecordStatus::Unusable(r)
        }
    }

    pub fn is_empty(&self) -> bool {
        self.partitions.is_empty()
    }

    pub fn markduplicates(&mut self) {
        //let mut dups = Vec::new();
        for (_f, tree) in self.partitions.drain() {
            //dups.clear();

            for tree_partition in tree.partitions(1) {
                let tree_index_most_occ = tree_partition[0];
                // no possible duplicates, one record in tree
                if tree[tree_index_most_occ].len() == 1 && tree_partition.len() == 1 {
                    let record_index = tree[tree_index_most_occ][0].index;
                    let name = self.current_reads[record_index].read_name().to_owned();
                    self.seen.insert(
                        name,
                        MarkResult {
                            corrected_umi: None,
                            is_dup: false,
                        },
                    );
                    continue;
                }

                let selected_umi = &tree[tree_index_most_occ][0].umi;
                let mut best = 0;
                let mut high_score = 0;

                let locs = tree_partition.iter()
                    .flat_map(|&ti| tree[ti].iter().map(|e| e.index))
                    .map(|i| self.current_reads[i].location.take().unwrap());
                let optical_dups = DuplicateClusters::new(locs);
                self.optical_dup_count += optical_dups.count_optical_dups(2500);
                eprintln!("cluster opt {}", self.optical_dup_count);


                for read_index in tree_partition
                    .into_iter()
                    .flat_map(|ti| tree[ti].iter().map(|e| e.index))
                {
                    let r = &mut self.current_reads[read_index];
                    let score = r.score();
                    if score > high_score {
                        best = read_index;
                        high_score = score;
                    }
                    let corrected_umi = if &r.umi != selected_umi {
                        r.correct_barcode_tag(String::from_utf8_lossy(selected_umi));
                        Some(selected_umi.to_vec())
                    } else {
                        None
                    };
                    r.flag_dup();
                    if r.is_paired() { 
                        self.seen.insert(
                            r.read_name().to_owned(),
                            MarkResult {
                                corrected_umi,
                                is_dup: true,
                            },
                        );
                    }
                }

                let r = &mut self.current_reads[best];
                r.unflag_dup();

                if r.is_paired() {
                    let res = self.seen.get_mut(r.read_name()).unwrap();
                    res.is_dup = false;
                }
            }
        }
    }

    pub fn clear_partitions(&mut self) {
        self.current_reads.clear();
        self.partitions.clear();
        self.in_partitions.clear();
        self.current_read_mates.clear();
        self.optical_dup_count = 0;
    }

    pub fn iter_records(&mut self) -> impl Iterator<Item = &UmiRecord> {
        for mate in &mut self.current_read_mates {
            let result = self.seen.remove(mate.read_name()).unwrap();
            if result.is_dup {
                mate.flag_dup();
            }
            if let Some(umi) = result.corrected_umi {
                mate.correct_barcode_tag(String::from_utf8_lossy(&umi));
            }
        }
        self.current_reads
            .iter()
            .chain(self.current_read_mates.iter())
    }
}
