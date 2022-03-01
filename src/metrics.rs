use std::fmt;

/// Duplication metrics.
#[derive(Debug, Default)]
pub struct Metrics {
    unpaired_reads_examined: usize,
    reads_pairs_examined: usize,
    secondary_or_supplementary_rds: usize,
    unmapped_reads: usize,
    unpaired_read_duplicates: usize,
    unpaired_read_optical_duplicates: usize,
    read_pair_duplicates: usize,
    read_pair_optical_duplicates: usize,
    corrected_umis: usize,
}

pub enum Status {
    UnpairedRead,
    ReadPair,
    SecondaryOrSupplementary,
    Unmapped,
    UnpairedDuplicate,
    UnpairedOpticalDuplicate,
    ReadpairDuplicate,
    ReadpairOpticalDuplicate,
    CorrectedUmi,
}

impl fmt::Display for Metrics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "UNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tCORRECTED_UMIS\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE")?;
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.unpaired_reads_examined,
            self.reads_pairs_examined,
            self.secondary_or_supplementary_rds,
            self.unmapped_reads,
            self.unpaired_read_duplicates,
            self.read_pair_duplicates,
            self.read_pair_optical_duplicates,
            self.corrected_umis,
            self.percent_duplication(),
            self.estimated_library_size()
        )
    }
}

impl Metrics {
    pub fn percent_duplication(&self) -> f32 {
        0.0
    }

    pub fn estimated_library_size(&self) -> u64 {
        0
    }

    pub fn count(&mut self,status : Status) {
        self.count_many(status, 1);
    }

    pub fn count_many(&mut self,status : Status, count: usize) {
        match status {
            Status::UnpairedRead => self.unpaired_reads_examined += count,
            Status::ReadPair => self.reads_pairs_examined += count,
            Status::SecondaryOrSupplementary => self.secondary_or_supplementary_rds += count,
            Status::Unmapped => self.unmapped_reads += count,
            Status::UnpairedDuplicate => self.unpaired_read_duplicates += count,
            Status::UnpairedOpticalDuplicate => self.unpaired_read_optical_duplicates += count,
            Status::ReadpairDuplicate => self.read_pair_duplicates += count,
            Status::ReadpairOpticalDuplicate => self.read_pair_optical_duplicates += count,
            Status::CorrectedUmi => self.corrected_umis += count,
        }
    }
}
