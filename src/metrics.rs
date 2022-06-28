use std::fmt;

use crate::record::Flags;

/// Duplication metrics.
#[derive(Debug, Default)]
pub struct Metrics {
    unpaired_reads_examined: usize,
    paired_reads_examined: usize,
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
    PairedRead,
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
        writeln!(f, "UNPAIRED_READS_EXAMINED\tPAIRED_READS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tUNPAIRED_READ_OPTICAL_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tCORRECTED_UMIS\tFRACTION_DUPLICATION\tESTIMATED_LIBRARY_SIZE")?;
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.unpaired_reads_examined,
            self.paired_reads_examined,
            self.secondary_or_supplementary_rds,
            self.unmapped_reads,
            self.unpaired_read_duplicates,
            self.unpaired_read_optical_duplicates,
            self.read_pair_duplicates,
            self.read_pair_optical_duplicates,
            self.corrected_umis,
            self.fraction_duplication(),
            self.estimated_library_size()
        )
    }
}

#[inline]
fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

impl Metrics {
    pub fn fraction_duplication(&self) -> f32 {
        (self.read_pair_duplicates / 2 + self.unpaired_read_duplicates) as f32
            / (self.unpaired_reads_examined + self.paired_reads_examined / 2) as f32
    }

    /// Library size estimation ported from picard DuplicationMetrics
    /// https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/markduplicates/EstimateLibraryComplexity.java
    ///
    /// modified to count single end reads as well
    pub fn estimated_library_size(&self) -> u64 {
        let read_pairs = ((self.paired_reads_examined - self.read_pair_optical_duplicates) / 2
            + self.unpaired_reads_examined
            - self.unpaired_read_optical_duplicates) as f64;
        let unique_read_pairs = ((self.paired_reads_examined - self.read_pair_duplicates) / 2
            + self.unpaired_reads_examined
            - self.unpaired_read_duplicates) as f64;
        let read_pair_duplicates = (read_pairs - unique_read_pairs) as f64;

        if read_pairs > 0.0 && read_pair_duplicates > 0.0 {
            let mut m = 1.0f64;
            let mut mm = 100.0f64;

            if unique_read_pairs >= read_pairs
                || f(m * unique_read_pairs, unique_read_pairs, read_pairs) < 0.0
            {
                panic!(
                    "Invalid values for pairs and unique pairs: {}, {}",
                    read_pairs, unique_read_pairs
                );
            }

            // find value of mm, large enough to act as other side for bisection method
            while f(mm * unique_read_pairs, unique_read_pairs, read_pairs) > 0.0 {
                mm *= 10.0;
            }

            // use bisection method (no more than 40 times) to find solution
            for _ in 0..40 {
                let r = (m + mm) / 2.0;
                let u = f(r * unique_read_pairs, unique_read_pairs, read_pairs);
                if u == 0.0 {
                    break;
                } else if u > 0.0 {
                    m = r;
                } else if u < 0.0 {
                    mm = r;
                }
            }

            (unique_read_pairs * (m + mm) / 2.0) as u64
        } else {
            0
        }
    }

    pub fn count_flags(&mut self, flags: Flags) {
        if flags.is_supplementary() || flags.is_secondary() {
            self.count(Status::SecondaryOrSupplementary);
        } else {
            if flags.is_unmapped() {
                self.count(Status::Unmapped);
            }

            if flags.is_segmented() {
                self.count(Status::PairedRead);
            } else {
                self.count(Status::UnpairedRead);
            }
        }
    }
    pub fn count_duplicate(&mut self, is_segmented: bool, is_optical: bool) {
        match (is_segmented, is_optical) {
            (true, true) => {
                self.count(Status::ReadpairOpticalDuplicate);
                self.count(Status::ReadpairDuplicate);
            }
            (true, false) => self.count(Status::ReadpairDuplicate),
            (false, true) => {
                self.count(Status::UnpairedDuplicate);
                self.count(Status::UnpairedOpticalDuplicate);
            }
            (false, false) => self.count(Status::UnpairedDuplicate),
        }
    }

    pub fn count(&mut self, status: Status) {
        self.count_many(status, 1);
    }

    pub fn count_many(&mut self, status: Status, count: usize) {
        match status {
            Status::UnpairedRead => self.unpaired_reads_examined += count,
            Status::PairedRead => self.paired_reads_examined += count,
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
