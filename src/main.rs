//! Filters records in a BAM that have a unique alignment and umi.

use anyhow::Result;

mod app;
pub mod bktree;
mod header;
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
