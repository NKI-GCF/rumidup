use std::convert::TryFrom;

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use noodles_bam::record::Record as NoodlesRecord;
use noodles_bam::{AsyncReader as NoodlesAsyncReader, AsyncWriter as NoodlesAsyncWriter};
use noodles_bgzf::AsyncReader as BgzfAsyncReader;
use noodles_bgzf::AsyncWriter as BgzfAsyncWriter;
use noodles_sam::header::{self, Header as BamHeader, ReferenceSequences};
use tokio::io::{AsyncRead, AsyncWrite};
use umibk::BkTree;

use super::record::{FragmentCoord, UmiRecord};
use crate::optical::*;

pub struct Reader<R>
where
    R: AsyncRead,
{
    bam: NoodlesReader<R>,
    record_buf: NoodlesRecord,
    header: String,
    pub reference_sequences: ReferenceSequences,
}

type NoodlesReader<R> = NoodlesAsyncReader<BgzfAsyncReader<R>>;
type NoodlesWriter<W> = NoodlesAsyncWriter<BgzfAsyncWriter<W>>;

impl<R> Reader<R>
where
    R: AsyncRead + std::marker::Unpin,
{
    pub async fn new(r: R) -> Result<Reader<R>> {
        let mut bam = NoodlesAsyncReader::new(r);
        let header = bam.read_header().await?;
        let reference_sequences = bam.read_reference_sequences().await?;

        Ok(Reader {
            bam,
            record_buf: NoodlesRecord::default(),
            header,
            reference_sequences,
        })
    }

    pub async fn read_record(&mut self) -> Result<Option<UmiRecord>> {
        let mut record = NoodlesRecord::default();
        match self.bam.read_record(&mut record).await? {
            0 => Ok(None),
            _n => Ok(Some(UmiRecord::try_from(record)?)),
        }
    }

    pub async fn out_bam<W: AsyncWrite + std::marker::Unpin>(
        &self,
        w: W,
    ) -> Result<NoodlesWriter<W>> {
        let mut writer = NoodlesAsyncWriter::new(w);
        let mut header: BamHeader = self.header.parse()?;
        header
            .programs_mut()
            .insert("umidedup".to_string(), header::Program::new("umidedup"));
        writer.write_header(&header).await?;
        writer
            .write_reference_sequences(&self.reference_sequences)
            .await?;

        Ok(writer)
    }
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

impl umibk::Dist for UmiIndex {
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
