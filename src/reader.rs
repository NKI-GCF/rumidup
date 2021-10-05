use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;

use anyhow::Result;
use itertools::Itertools;
use noodles_bam::{AsyncReader as NoodlesAsyncReader, AsyncWriter as NoodlesAsyncWriter};
use noodles_bam::record::Record as NoodlesRecord;
use noodles_sam::header::{self, Header as BamHeader, ReferenceSequences};
use noodles_sam::record::data::field::tag::Tag;
use noodles_sam::record::{Flags, data::field::Field, data::field::value::Value};
use tokio::io::{AsyncRead, AsyncWrite};

use super::record::{UmiRecord, FragmentCoord};
use umibk::BkTree;

pub struct Reader<R> where R: AsyncRead {
    bam: NoodlesAsyncReader<R>,
    record_buf: NoodlesRecord,
    header: String,
    pub reference_sequences: ReferenceSequences, 
}

impl<R> Reader<R> where R: AsyncRead + std::marker::Unpin {
    pub async fn new(r: R) -> Result<Reader<R>> {
        let mut bam = NoodlesAsyncReader::new(r);
        let header = bam.read_header().await?;
        let reference_sequences = bam.read_reference_sequences().await?;

        Ok(Reader { bam, record_buf: NoodlesRecord::default(), header, reference_sequences })
    }

    pub async fn read_record(&mut self) -> Result<Option<UmiRecord>> {
        match self.bam.read_record(&mut self.record_buf).await? {
            0 => Ok(None),
            _n => {
                let record = self.record_buf.try_into_sam_record(&self.reference_sequences)?;
                Ok(Some(UmiRecord::try_from(record)?))
            }
        }
    }

    pub async fn out_bam<W: AsyncWrite + std::marker::Unpin>(&self, w: W) -> Result<NoodlesAsyncWriter<W>> {
        let mut writer = NoodlesAsyncWriter::new(w);
        let mut header: BamHeader = self.header.parse()?;
        header.programs_mut().insert("umidedup".to_string (), header::Program::new("umidedup"));
        writer.write_header(&header).await?;
        writer.write_reference_sequences(&self.reference_sequences).await?;

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

struct Partition {
    records: BkTree<UmiIndex>,
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
            self.umi.iter().zip(b.umi.iter()).filter(|&(a, b)| a == &b'N' || b == &b'N' || a != b).count()
        } else {
            std::cmp::max(la, lb)
        }
    }
}

pub struct ReadPartitions {
    pub seen: HashMap<String, bool>,
    pub current_reads: Vec<UmiRecord>,
    partitions: HashMap<(FragmentCoord, FragmentCoord), BkTree<UmiIndex>>,
}

pub enum RecordStatus {
    Unusable(UmiRecord),
    PossibleDup,
    Duplicate,
    SeenMate(UmiRecord)
}

impl ReadPartitions {
    pub fn new() -> ReadPartitions {
        ReadPartitions { seen: HashMap::new(), current_reads: Vec::new(), partitions: HashMap::new() }
    }

    pub fn add_record(&mut self, mut r: UmiRecord) -> RecordStatus {
        if r.is_paired() {
            if let Some(&is_dup) = self.seen.get(r.read_name()) {
                if is_dup {
                    r.flag_dup();
                }
                self.seen.remove(r.read_name());
                return RecordStatus::SeenMate(r);
            }
        }

        //determine fragment coordinates
        if let Some(fragment) = r.fragment_markers() {
            let record_idx = self.current_reads.len();
            let umi = UmiIndex { umi: r.umi.to_vec(), index: record_idx};
            self.current_reads.push(r);

            //insert into tree
            let aap = fragment.clone();
            let e = self.partitions.entry(fragment).or_insert_with(BkTree::new);
            if e.insert(umi) {
            if self.current_reads.last().unwrap().record.position().map(|p| p.into()) == Some(52868755) {
                eprintln!("{:?} {:?}", aap, self.partitions.get(&aap));
            }
                RecordStatus::PossibleDup
            } else {
            if self.current_reads.last().unwrap().record.position().map(|p| p.into()) == Some(52868755) {
                eprintln!("{:?} {:?}", aap, self.partitions.get(&aap));
            }
                RecordStatus::Duplicate
            }
        } else {
            return RecordStatus::Unusable(r);
        }

    }

    pub fn is_empty(&self) -> bool {
        self.partitions.is_empty()
    }

    pub fn markduplicates(&mut self, d: bool) {

        //let mut dups = Vec::new();
        for (_f, tree) in self.partitions.drain() {
            //dups.clear();

            for tree_partition in tree.partitions(1) {
                let tree_index_most_occ = tree_partition[0];
                if tree[tree_index_most_occ].len() == 1 && tree_partition.len() == 1 {
                    continue
                }
                
                /*
                let indexes: Vec<_> = tree_partition.into_iter()
                    .flat_map(|ti| tree[ti].iter().map(|e| e.index))
                    .collect();
*/
                let corrected_umi =  &tree[tree_index_most_occ][0].umi;
                let mut best = 0;
                let mut high_score = 0;

                for read_index in  tree_partition.into_iter().flat_map(|ti| tree[ti].iter().map(|e| e.index)) {
                    let r = &mut self.current_reads[read_index];
                    let score = r.score();
                    if score > high_score {
                        best = read_index;
                        high_score = score;
                    }
                    if &r.umi != corrected_umi {
                        r.modify_barcode_tag(String::from_utf8_lossy(&corrected_umi));
                    }
                    r.flag_dup();
                    if r.is_paired() {
                        self.seen.insert(r.read_name().to_owned(), true);
                    }
                }

                let r = &mut self.current_reads[best];
                r.unflag_dup();
               
                if r.is_paired() {
                    *self.seen.get_mut(r.read_name()).unwrap() = false;
                }
            }
        }
    }

    pub fn clear_partitions(&mut self) {
        self.current_reads.clear();
        self.partitions.clear();
    }
}

fn processdups(record_indices: &[usize], records: &mut [UmiRecord]) {
    let final_umi = records[record_indices[0]].umi.to_vec();

    for &i in record_indices.iter().skip(1) {
        let r = &mut records[i];
        if r.umi != final_umi {
            r.modify_barcode_tag(String::from_utf8_lossy(&final_umi));
        }
        r.flag_dup();
    }
}
