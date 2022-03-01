//! Filters records in a BAM that have a unique alignment and umi.
//!

//e.g. /net/NGSanalysis/data/m.hoekstra/6434/6434_8_2h_1_ng_ml_IFNg_10_ng_ml_TNFa_CCAGTCGTCA-CGCTAGGCTA_S25.bam

use anyhow::Result;
use tokio;

pub mod io;
pub mod record;
pub mod optical;
pub mod bktree;

use io::RecordStatus;

#[tokio::main(flavor = "current_thread")]
async fn main() -> Result<()> {
    let stdin = tokio::io::stdin();
    let stdout = tokio::io::stdout();

    let mut io = io::BamIo::new(stdin, stdout).await?;

    //let mut pos = None;
    let mut parts = io::ReadPartitions::default();
    let mut optical_dups = 0;
    /*
    while let Some(mut r) = io.read_record().await? {
        //unmark current setting
        r.unflag_dup();
        //steps
        // is read mapped?
        // already marked?
        // is mate mapped?
        // is paired?
        // is primary alignment
        //FIXME also check read chr changes
        let read_pos = r.position();
        if read_pos.is_some() && read_pos != pos {
            if !parts.is_empty() {
                //dedups parts
                if parts.current_reads.len() > 50 {
                    eprintln!("dedup bundle for pos {:?}", read_pos);
                }
                parts.markduplicates();
                optical_dups += parts.optical_dup_count;
                for r in parts.iter_records() {
                    writer.write_record(&r.record).await?;
                }
                //dump mates using their mate result
                parts.clear_partitions();
            }
            pos = read_pos;
        }

        match parts.add_record(r) {
            RecordStatus::PossibleDup | RecordStatus::Duplicate | RecordStatus::InPartition => {}
            RecordStatus::Unusable(r) | RecordStatus::SeenMate(r) => {
                reader.write_record(&r.into()).await?;
            }
        }
    }

        */
    let mut bundle = Vec::new();
    while io.read_bundle(&mut bundle).await? {

        if bundle.len() == 1 {
            for r in bundle.drain(..) {
                io.write_record(&r).await?;
            }
        } else {
            if bundle.len() > 50 {
                eprintln!("Big bundle ({}) at {:?} {:?}", bundle.len(), bundle[0].reference_sequence_id(), bundle[0].position());
            }
            for r in bundle.drain(..) {
                let r = record::UmiRecord::try_from(r).unwrap();
                match parts.add_record(r) {
                    RecordStatus::PossibleDup | RecordStatus::Duplicate | RecordStatus::InPartition => {}
                    RecordStatus::Unusable(r) | RecordStatus::SeenMate(r) => {
                        //writer.write_record(&r.into()).await?;
                    }
                }
            }
            parts.markduplicates();
            optical_dups += parts.optical_dup_count;
            for r in parts.iter_records() {
                io.write_record(&r.record).await?;
            }
            parts.clear_partitions();
        }
    }

    if !parts.seen.is_empty() {
        eprintln!("{} mates not found", parts.seen.len());
    }

    eprintln!("opt dups {}", optical_dups);

    io.shutdown().await?;

    Ok(())
}
