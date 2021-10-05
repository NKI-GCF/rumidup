//! Filters records in a BAM that have a unique alignment and umi.
//!

//e.g. /net/NGSanalysis/data/m.hoekstra/6434/6434_8_2h_1_ng_ml_IFNg_10_ng_ml_TNFa_CCAGTCGTCA-CGCTAGGCTA_S25.bam

use anyhow::Result;
use tokio::io;
use clap::App;


pub mod record;
pub mod reader;

use reader::RecordStatus;

#[tokio::main(flavor="current_thread")]
async fn main() -> Result<()> {

    let stdin = io::stdin();
    let stdout = io::stdout();

    let mut reader = reader::Reader::new(stdin).await?;
    let mut writer = reader.out_bam(stdout).await?;

    let mut pos = None;
    let mut parts = reader::ReadPartitions::new();
    while let Some(r) = reader.read_record().await? {
        //steps
        // is read mapped?
        // already marked?
        // is mate mapped?
        // is paired?
        // is primary alignment
        let read_pos = r.position();
        if read_pos.is_some() && read_pos != pos && !parts.is_empty() {
            //dedups parts
            if parts.current_reads.len() > 50 {
                eprintln!("dedup bundle for pos {:?}", read_pos);
            }
            parts.markduplicates(parts.current_reads.len() > 50);
            for r in parts.current_reads.drain(..) {
                writer.write_sam_record(&reader.reference_sequences, &r.record.into()).await?;
            }
            pos = read_pos;
            parts.clear_partitions();
        }

        match parts.add_record(r) {
            RecordStatus::PossibleDup | RecordStatus::Duplicate => {},
            RecordStatus::Unusable(r) | RecordStatus::SeenMate(r) => {
                writer.write_sam_record(&reader.reference_sequences, &r.into()).await?;
            }
        }
    }

    if !parts.seen.is_empty() {
        eprintln!("{} mates not found", parts.seen.len());
    }

    writer.shutdown().await?;

    Ok(())

}

