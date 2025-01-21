use anyhow::Result;
use noodles_bam::{AsyncReader as NoodlesAsyncReader, AsyncWriter as NoodlesAsyncWriter};
use noodles_bgzf::{
    writer::CompressionLevel, AsyncReader as BgzfAsyncReader, AsyncWriter as BgzfAsyncWriter,
};
use noodles_core::position::Position;
use noodles_sam::alignment::record_buf::RecordBuf as NoodlesRecord;
use thiserror::Error;
use tokio::io::{self, AsyncRead, AsyncWrite};

use crate::app::CmdInfo;
use crate::header::{self, BamHeader};

pub struct BamIo<R, W>
where
    R: AsyncRead,
    W: AsyncWrite,
{
    in_bam: NoodlesReader<R>,
    out_bam: NoodlesWriter<W>,
    next_record: Option<NoodlesRecord>,
    header: BamHeader,
    //pub reference_sequences: ReferenceSequences,
}

type NoodlesReader<R> = NoodlesAsyncReader<BgzfAsyncReader<R>>;
type NoodlesWriter<W> = NoodlesAsyncWriter<BgzfAsyncWriter<W>>;

// Combine reference id and position to group reads
type ChrPos = (usize, Position);

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
    /// Pass a CmdInfo struct to insert PG line in the BAM header
    pub async fn new(
        read: R,
        write: W,
        force_rerun: bool,
        compress_out: bool,
        exe_info: Option<CmdInfo>,
    ) -> Result<BamIo<R, W>, BamIoError> {
        let mut in_bam = NoodlesAsyncReader::new(read);
        let mut header: BamHeader = in_bam.read_header().await?.into();
        //let reference_sequences = in_bam.read_reference_sequences().await?;

        //let mut header: BamHeader = header_string.parse()?;
        if let Some(pg) = header.detect_markdups() {
            eprintln!(
                "Probable previous markduplictes detected in PG lines:\n{}",
                pg
            );
            if !force_rerun {
                return Err(BamIoError::MarkDupsDetected);
            }
        }

        if let Some(cmd) = exe_info {
            header.add_rumidup_pg(&cmd.command_line, &cmd.version)
        }
        let mut builder = noodles_bgzf::r#async::writer::Builder::default();
        if !compress_out {
            builder = builder.set_compression_level(CompressionLevel::none());
        }
        let writer = builder.build_with_writer(write);
        let mut out_bam = NoodlesAsyncWriter::from(writer);
        out_bam.write_header(header.as_ref()).await?;
        /*
        out_bam
            .write_reference_sequences(&reference_sequences)
            .await?;
        */
        Ok(BamIo {
            in_bam,
            out_bam,
            next_record: None,
            header,
            //reference_sequences,
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
    pub async fn read_at_pos(
        &mut self,
        chr_pos: Option<ChrPos>,
    ) -> io::Result<Option<NoodlesRecord>> {
        if self.next_record.is_none() {
            if let Some(record) = self.read_from_bam().await? {
                self.next_record = Some(record);
            } else {
                return Ok(None);
            }
        }

        let next_chr_pos = self.get_chr_pos(self.next_record.as_ref().unwrap());
        if next_chr_pos == chr_pos {
            Ok(Some(self.next_record.take().unwrap()))
        } else {
            Ok(None)
        }
    }

    /// Read a new record from the stream and return it.
    async fn read_from_bam(&mut self) -> io::Result<Option<NoodlesRecord>> {
        let mut record = NoodlesRecord::default();
        match self
            .in_bam
            .read_record_buf(self.header.as_ref(), &mut record)
            .await?
        {
            0 => Ok(None),
            _n => Ok(Some(record)),
        }
    }
    fn get_chr_pos(&self, record: &NoodlesRecord) -> Option<ChrPos> {
        record.reference_sequence_id().zip(record.alignment_start())
    }

    pub async fn read_bundle(&mut self, bundle: &mut Vec<NoodlesRecord>) -> io::Result<bool> {
        bundle.clear();
        if let Some(first) = self.read().await? {
            let chr_pos = self.get_chr_pos(&first);
            bundle.push(first);
            if chr_pos.is_some() {
                while let Some(next) = self.read_at_pos(chr_pos).await? {
                    bundle.push(next);
                }
            }
        }
        Ok(!bundle.is_empty())
    }

    pub async fn read_bundle_owned(&mut self) -> io::Result<Vec<NoodlesRecord>> {
        let mut bundle = Vec::new();
        if let Some(first) = self.read().await? {
            let chr_pos = self.get_chr_pos(&first);
            bundle.push(first);
            if chr_pos.is_some() {
                while let Some(next) = self.read_at_pos(chr_pos).await? {
                    bundle.push(next);
                }
            }
        }
        Ok(bundle)
    }

    pub async fn write_record(&mut self, record: &NoodlesRecord) -> io::Result<()> {
        self.out_bam
            .write_alignment_record(self.header.as_ref(), record)
            .await
    }

    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.out_bam.shutdown().await
    }
}

#[derive(Debug, Error)]
pub enum BamIoError {
    #[error("Io error reading BAM")]
    IoError(#[from] std::io::Error),
    #[error("ParseError reading BAM: {0}")]
    ParseError(String),
    #[error("Error in BAM header")]
    HeaderError(#[from] header::ParseError),
    #[error("Markduplicates evidence present in BAM header")]
    MarkDupsDetected,
}
