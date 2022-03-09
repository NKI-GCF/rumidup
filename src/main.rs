//! Filters records in a BAM that have a unique alignment and umi.
//!

//e.g. /net/NGSanalysis/data/m.hoekstra/6434/6434_8_2h_1_ng_ml_IFNg_10_ng_ml_TNFa_CCAGTCGTCA-CGCTAGGCTA_S25.bam

use anyhow::Result;

mod app;
pub mod bktree;
pub mod io;
pub mod markdups;
pub mod metrics;
pub mod optical;
pub mod record;

use app::App;

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let mut app = App::new().await?;
    app.run().await?;

    Ok(())
}
