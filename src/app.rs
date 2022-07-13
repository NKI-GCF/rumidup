use std::marker::Unpin;
use std::path::PathBuf;

use ahash::{AHashMap, AHashSet};
use clap::Parser;
use thiserror::Error;
use tokio::{
    fs::File,
    io::{self, AsyncRead, AsyncWrite},
};

use crate::{
    io::{BamIo, BamIoError},
    markdups::{CorrectedUmi, Dist, MarkResult, UmiClusters},
    metrics::{Metrics, Status},
    optical::DuplicateClusters,
    record::{BamRecord, ReadName, RecordError, UmiRecord},
};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Config {
    /// The input bam file. rumidup reads from stdin when omitted
    #[clap(short, long)]
    pub bam: Option<PathBuf>,

    /// The output bam file. rumidup writes to stdout when omitted
    #[clap(short, long)]
    pub output: Option<PathBuf>,

    /// The duplication metrics file, if missing metrics will be written to stderr
    #[clap(short = 'm', long)]
    pub metrics: Option<PathBuf>,

    /// Ignore previous duplicate marking applied to BAM file. This information is extracted from
    /// the header. Use --force to redoing duplicate marking
    #[clap(short, long)]
    pub force: bool,

    /// When a UMI is corrected the original tag is writtento the OX field.
    /// use --no-original-tag to suppress this
    #[clap(short = 'x', long)]
    pub no_original_tag: bool,

    /// Don't remove UMI from read name
    #[clap(short = 'k', long)]
    pub keep_readname: bool,

    /// UMI distance. The maximin hamming distance between the UMI seqeunces used to
    /// consider read(pairs) to be duplicates
    #[clap(short = 'd', long, default_value = "1")]
    pub umi_distance: usize,

    /// Optical duplicate pixel distance.
    /// Maximum distance between clusters to consider them optical duplicates
    /// use 100 for HiSeq/NextSeq, 2500 for NovaSeq.
    /// use 0 to disable optical duplicate counting
    /// Only affects metrics
    #[clap(short = 'p', long, default_value = "100")]
    pub pixel_distance: i32,
}

pub struct App {
    config: Config,
    bamio: BamIo<Box<dyn AsyncRead + Unpin>, Box<dyn AsyncWrite + Unpin>>,
    seen: AHashMap<ReadName, MarkResult>,
    metrics: Metrics,
    records: Vec<BamRecord>,
    umirecords: Vec<UmiRecord>,
    record_results: Vec<MarkResult>,
}

/// In the BK tree we store the umi, and the index of the
/// record in the ReadBundle
#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd)]
struct UmiIndex {
    index: usize,
    umi: Vec<u8>,
}

/// The Dist trait impl for the BkTree
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

impl App {
    pub async fn new() -> Result<App, RumidupError> {
        let config = Config::parse();

        let read: Box<dyn AsyncRead + Unpin> = if let Some(p) = config.bam.as_ref() {
            Box::new(File::open(p).await?)
        } else {
            Box::new(io::stdin())
        };

        let write: Box<dyn AsyncWrite + Unpin> = if let Some(p) = config.output.as_ref() {
            Box::new(File::create(p).await?)
        } else {
            Box::new(io::stdout())
        };

        let bamio = BamIo::new(read, write).await?;

        Ok(App {
            config,
            bamio,
            seen: AHashMap::new(),
            metrics: Metrics::default(),
            records: Vec::new(),
            umirecords: Vec::new(),
            record_results: Vec::new(),
        })
    }

    /// Deduplication is performed bundles of records that all map to the same start location.
    /// For paired-end reads only one read is used to determine the duplcation status.
    pub async fn run(&mut self) -> Result<(), RumidupError> {
        while self.read_bundle().await? {
            self.process_in_records()?;
            self.mark_unknown()?;
            self.update_records()?;
            self.write_records().await?;
        }

        self.bamio.shutdown().await?;
        self.write_report()?;

        Ok(())
    }

    async fn read_bundle(&mut self) -> Result<bool, std::io::Error> {
        self.bamio.read_bundle(&mut self.records).await
    }

    fn process_in_records(&mut self) -> Result<(), RumidupError> {
        for record in self.records.drain(..) {
            let mut umirecord: UmiRecord = record.into();
            let mut result = MarkResult::Unusable;

            let flags = umirecord.flags();
            self.metrics.count_flags(flags);

            umirecord.extract_umi(!self.config.keep_readname)?;

            //duplication status of reads that have both (first and last segment) mapped is only
            //determined once and applied to both
            let is_mapped_pair = !flags.is_secondary()
                && !flags.is_supplementary()
                && flags.is_segmented()
                && !flags.is_unmapped()
                && !flags.is_mate_unmapped();

            //check if mate read has been processed before
            let seen = is_mapped_pair && self.seen.contains_key(umirecord.read_name());

            let willmarkdup =
                !seen && !flags.is_secondary() && !flags.is_supplementary() && !flags.is_unmapped();

            if willmarkdup && is_mapped_pair {
                umirecord.extract_mate_tags()?;
            }

            if willmarkdup && self.config.pixel_distance > 0 {
                umirecord.extract_location()?;
            }

            if seen {
                result = MarkResult::Skipped;
            } else if willmarkdup {
                result = if let Some(coord) = umirecord.fragment_markers() {
                    MarkResult::Unknown(coord)
                } else {
                    MarkResult::Failed
                };

                if is_mapped_pair {
                    self.seen
                        .insert(umirecord.read_name().clone(), result.clone());
                }
            }
            self.record_results.push(result);
            self.umirecords.push(umirecord);
        }

        Ok(())
    }

    fn mark_unknown(&mut self) -> Result<(), RumidupError> {
        //extract status unknown records and create umi trees
        let umi_clusters: UmiClusters<_> = self
            .umirecords
            .iter_mut()
            .enumerate()
            .zip(self.record_results.iter())
            .filter_map(|((index, r), s)| match s {
                MarkResult::Unknown(coord) => Some((
                    *coord,
                    UmiIndex {
                        index,
                        umi: r.umi.clone(),
                    },
                )),
                _ => None,
            })
            .collect();

        for coord in umi_clusters.coordinate_pairs() {
            for dups in umi_clusters.mark_duplicates(coord, self.config.umi_distance) {
                if dups.len() == 1 {
                    self.set_record_result(
                        dups[0].index,
                        MarkResult::Unique(CorrectedUmi::default()),
                    );
                    continue;
                }

                //extract the read with the highers baseq
                let best = (0..dups.len())
                    .max_by_key(|&i| self.umirecords[dups[i].index].score())
                    .unwrap();

                //mark the best record as unique
                self.set_record_result(
                    dups[best].index,
                    MarkResult::Unique(CorrectedUmi::from_cmp(&dups[best].umi, &dups[0].umi)),
                );

                if self.config.pixel_distance == 0 {
                    for uix in dups.iter().filter(|i| i.index != best) {
                        self.set_record_result(
                            uix.index,
                            MarkResult::Duplicate(CorrectedUmi::from_cmp(&uix.umi, &dups[0].umi)),
                        );
                    }
                } else {
                    //count opticals
                    let opticals = DuplicateClusters::new(
                        dups.iter()
                            .map(|ui| self.umirecords[ui.index].location.take().unwrap()),
                    );
                    let optclust_ids = opticals.optical_cluster_ids(self.config.pixel_distance);
                    let mut cluster_ids = AHashSet::new();
                    cluster_ids.insert(optclust_ids[best]);

                    for (index, cluster) in optclust_ids
                        .into_iter()
                        .enumerate()
                        .filter(|&(index, _optid)| index != best)
                    {
                        if cluster_ids.insert(cluster) {
                            //mark as pcr dup
                            self.set_record_result(
                                dups[index].index,
                                MarkResult::Duplicate(CorrectedUmi::from_cmp(
                                    &dups[index].umi,
                                    &dups[0].umi,
                                )),
                            );
                        } else {
                            // mark as optical dup
                            self.set_record_result(
                                dups[index].index,
                                MarkResult::OpticalDuplicate(CorrectedUmi::from_cmp(
                                    &dups[index].umi,
                                    &dups[0].umi,
                                )),
                            );
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn set_record_result(&mut self, i: usize, result: MarkResult) {
        if self.umirecords[i].is_paired() {
            self.seen
                .insert(self.umirecords[i].read_name().to_owned(), result.clone());
        }
        self.record_results[i] = result;
    }

    fn update_records(&mut self) -> Result<(), RumidupError> {
        let mut old_result;
        for (i, mut result) in self.record_results.iter().enumerate() {
            if matches!(result, MarkResult::Skipped) {
                old_result = self.seen.remove(self.umirecords[i].read_name());
                if old_result.is_some() {
                    result = old_result.as_ref().unwrap();
                } else {
                    eprintln!("No old for mate seen {:?}", self.seen);
                }
            }
            match result {
                MarkResult::Duplicate(umi) => {
                    if umi.is_corrected() {
                        self.umirecords[i]
                            .correct_barcode_umi(umi.to_string(), !self.config.no_original_tag);
                        self.metrics.count(Status::CorrectedUmi);
                    }
                    self.umirecords[i].flag_dup();
                    self.metrics
                        .count_duplicate(self.umirecords[i].is_paired(), false);
                }
                MarkResult::OpticalDuplicate(umi) => {
                    if umi.is_corrected() {
                        self.umirecords[i]
                            .correct_barcode_umi(umi.to_string(), !self.config.no_original_tag);
                        self.metrics.count(Status::CorrectedUmi);
                    }
                    self.umirecords[i].flag_dup();
                    self.metrics
                        .count_duplicate(self.umirecords[i].is_paired(), true);
                }
                MarkResult::Unique(umi) => {
                    if umi.is_corrected() {
                        self.umirecords[i]
                            .correct_barcode_umi(umi.to_string(), !self.config.no_original_tag);
                        self.metrics.count(Status::CorrectedUmi);
                    }
                    self.umirecords[i].unflag_dup();
                }
                MarkResult::Unusable => {}
                MarkResult::Skipped => {
                    panic!("Nested skipped/unusable record {:?}", self.umirecords[i])
                }
                MarkResult::Failed => eprintln!(
                    "Record failed to parse usable coordinates {}, rumidup skipped this record, but it might indicate problems with the BAM file",
                    self.umirecords[i].read_name()
                ),
                MarkResult::Unknown(_) => eprintln!(
                    "Record still marked unknown after dedup: {}. This is a BUG!",
                    self.umirecords[i].read_name()
                ),
            }
        }

        Ok(())
    }

    /// empties the umirecords and the results Vec
    async fn write_records(&mut self) -> io::Result<()> {
        for record in self.umirecords.drain(..) {
            self.bamio.write_record(&record.into()).await?;
        }
        self.record_results.clear();
        Ok(())
    }

    fn write_report(&self) -> io::Result<()> {
        if let Some(path) = &self.config.metrics {
            use std::io::Write;
            let mut mout = std::fs::File::create(path)?;
            write!(mout, "{}", self.metrics)?;
        } else {
            eprintln!("{}", self.metrics);
        }

        Ok(())
    }
}

#[derive(Debug, Error)]
#[allow(clippy::enum_variant_names)]
pub enum RumidupError {
    #[error("IoError")]
    IoError(#[from] std::io::Error),
    #[error("Error r BAM record")]
    ReaderError(#[from] BamIoError),
    #[error("Error using BAM record")]
    RecordError(#[from] RecordError),
}
