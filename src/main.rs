//! Filters records in a BAM that have a unique alignment and umi.
//!

//e.g. /net/NGSanalysis/data/m.hoekstra/6434/6434_8_2h_1_ng_ml_IFNg_10_ng_ml_TNFa_CCAGTCGTCA-CGCTAGGCTA_S25.bam

use anyhow::Result;
use tokio::io;

pub mod reader;
pub mod record;
pub mod optical;

use reader::RecordStatus;

#[tokio::main(flavor = "current_thread")]
async fn main() -> Result<()> {
    let stdin = io::stdin();
    let stdout = io::stdout();

    let mut reader = reader::Reader::new(stdin).await?;
    let mut writer = reader.out_bam(stdout).await?;

    let mut pos = None;
    let mut parts = reader::ReadPartitions::default();
    let mut optical_dups = 0;
    while let Some(mut r) = reader.read_record().await? {
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
                writer.write_record(&r.into()).await?;
            }
        }
    }

    if !parts.seen.is_empty() {
        eprintln!("{} mates not found", parts.seen.len());
    }

    eprintln!("opt dups {}", optical_dups);

    writer.shutdown().await?;

    Ok(())
}
