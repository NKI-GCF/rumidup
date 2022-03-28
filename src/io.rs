use anyhow::Result;
use noodles_bam::record::Record as NoodlesRecord;
use noodles_bam::{AsyncReader as NoodlesAsyncReader, AsyncWriter as NoodlesAsyncWriter};
use noodles_bgzf::AsyncReader as BgzfAsyncReader;
use noodles_bgzf::AsyncWriter as BgzfAsyncWriter;
use noodles_sam::header::{self, Header as BamHeader, ReferenceSequences};
use thiserror::Error;
use tokio::io::{self, AsyncRead, AsyncWrite};

pub struct BamIo<R, W>
where
    R: AsyncRead,
    W: AsyncWrite,
{
    in_bam: NoodlesReader<R>,
    out_bam: NoodlesWriter<W>,
    next_record: Option<NoodlesRecord>,
    pub reference_sequences: ReferenceSequences,
    ref_shift: u32,
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

        let out_bam = out_bam(write, &header, &reference_sequences).await?;
        let ref_shift = reference_sequences
            .len()
            .next_power_of_two()
            .leading_zeros();

        Ok(BamIo {
            in_bam,
            out_bam,
            next_record: None,
            reference_sequences,
            ref_shift,
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
        chr_pos: Option<usize>,
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
        match self.in_bam.read_record(&mut record).await? {
            0 => Ok(None),
            _n => Ok(Some(record)),
        }
    }
    fn get_chr_pos(&self, record: &NoodlesRecord) -> Option<usize> {
        record.position().and_then(|p| {
            usize::try_from(i32::from(p))
                .expect("Negative position?")
                .checked_shl(self.ref_shift)
                .and_then(|p| record.reference_sequence_id().map(|c| usize::from(c) | p))
        })
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
        self.out_bam.write_record(record).await
    }

    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.out_bam.shutdown().await
    }
}

async fn out_bam<W: AsyncWrite + std::marker::Unpin>(
    w: W,
    header: &str,
    reference_sequences: &ReferenceSequences,
) -> Result<NoodlesWriter<W>, BamIoError> {
    let mut writer = NoodlesAsyncWriter::new(w);
    let mut header: BamHeader = header
        .parse()
        .map_err(|e: header::ParseError| BamIoError::ParseError(e.to_string()))?;
    header
        .programs_mut()
        .insert("rumidup".to_string(), header::Program::new("rumidup"));
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
}
