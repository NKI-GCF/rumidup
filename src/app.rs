use std::path::PathBuf;
use std::marker::Unpin;

use ahash::AHashMap;
use clap::Parser;
use noodles_sam as sam;
use thiserror::Error;
use tokio::{
    io::{self, AsyncRead, AsyncWrite},
    fs::File,
};

use crate::{
    io::{BamIo, BamIoError},
    markdups::{MarkResult},
    metrics::{Metrics, Status},
    record::{BamRecord, ReadName, RecordError},
};


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Config {
    /// The input bam file. umidedup reads from stdin when omitted
    #[clap(short, long)]
    pub bam: Option<PathBuf>,

    /// The output bam file. umidedup writes to stdout when omitted
    #[clap(short, long)]
    pub output: Option<PathBuf>,

    /// The duplication metrics file, if missing metrics will be written to stderr
    #[clap(short = 'm', long)]
    pub metrics: Option<PathBuf>,

    /// Don't touch reads already flagged as duplicate in input bam.
    /// By default every read(pair) is evaluated.
    #[clap(short = 's', long)]
    pub skip_marked_duplcates: bool,

    /// By default umidedup will try to parse the UMI from the readname or the RX tag
    /// You can specify a custom tag if you want to override this behaviour
    #[clap(short = 'u', long)]
    pub umi_tag: Option<sam::record::data::field::Tag>,

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
    pub pixel_distance: u32,
}

pub struct App {
    config: Config,
    bamio: BamIo<Box<dyn AsyncRead + Unpin>, Box<dyn AsyncWrite + Unpin>>,
    seen: AHashMap<ReadName, MarkResult>,
    metrics: Metrics,
}


impl App {
    pub async fn new() -> Result<App, UmiDedupError> {
        let config = Config::parse();

        let read:Box<dyn AsyncRead + Unpin> = if let Some(p) = config.bam.as_ref() {
            Box::new(File::open(p).await?)
        } else {
            Box::new(io::stdin())
        };

        let write:Box<dyn AsyncWrite + Unpin> = if let Some(p) = config.output.as_ref() {
            Box::new(File::create(p).await?)
        } else {
            Box::new(io::stdout())
        };

        let bamio = BamIo::new(read, write).await?;

        Ok(App { config, bamio, seen: AHashMap::new(), metrics: Metrics::default() })
    }

    /// Deduplication is performed bundles of records that all map to the same start location.
    /// For paired-end reads only one read is used to determine the duplcation status.
    pub async fn run(&mut self) -> Result<(), UmiDedupError> {

        let mut bundle = Vec::new();
        while self.bamio.read_bundle(&mut bundle).await? {
            for record in &bundle {
            }
            //extract records for deduplication
            // run dedup
            //
            // update metrics
            //
            // update records
            //
            // write bundle to out
        }

        // print/write the metrics (sync)
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
pub enum UmiDedupError {
    #[error("IoError")]
    IoError(#[from] std::io::Error),
    #[error("Error r BAM record")]
    ReaderError(#[from] BamIoError),
    #[error("Error using BAM record")]
    RecordError(#[from] RecordError),
}

