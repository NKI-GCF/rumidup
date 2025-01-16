use std::str::FromStr;

use ahash::AHashMap;
pub use noodles_sam::header::ParseError;
use noodles_sam::header::{Header, record::value::map::{Map, Program}};

pub struct BamHeader(Header);

impl From<Header> for BamHeader {
    fn from(h: Header) -> BamHeader {
        BamHeader(h)
    }
}

impl FromStr for BamHeader {
    type Err = ParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BamHeader(Header::from_str(s)?))
    }
}

impl AsRef<Header> for BamHeader {
    fn as_ref(&self) -> &Header {
        &self.0
    }
}

impl BamHeader {
    pub fn detect_markdups(&self) -> Option<String> {
        self.0
            .programs()
            .values()
            .find(|p| p.is_markdup())
            .map(|p| p.to_string())
    }

    pub fn add_rumidup_pg(&mut self, command_line: &str, version: &str) {
        let pps = self.pg_chains();

        if pps.is_empty() {
            let id = self.next_rumidup_id();
            self.0.programs_mut().insert(
                id.clone(),
                Map::<Program>::builder()
                    .set_name("rumidup")
                    .set_version(version)
                    .set_command_line(command_line)
                    .build()
                    .expect("Error building program"),
            );
        } else {
            for p in pps {
                let id = self.next_rumidup_id();
                self.0.programs_mut().insert(
                    id.clone(),
                    Map::<Program>::builder()
                        .set_name("rumidup")
                        .set_version(version)
                        .set_command_line(command_line)
                        .set_previous_id(p)
                        .build()
                        .expect("Error building program"),
                );
            }
        }
    }

    fn pg_chains(&self) -> Vec<String> {
        let programs = self.0.programs();

        let mut chains = AHashMap::new();
        let mut ids: Vec<&String> = programs.keys().collect();

        // extract all the chain heads (non or invalid PP)
        ids.retain(|&id| {
            if let Some(previous_id) = programs.get(id).unwrap().previous_id() {
                //check invalid PP ref
                if !programs.contains_key(previous_id) {
                    eprintln!(
                        "Error in PG chain. Previous id does not exist: {}",
                        previous_id
                    );
                    chains.insert(previous_id, 1);
                    false
                } else {
                    true
                }
            } else {
                chains.insert(id, 1);
                false
            }
        });

        // append remainder to head until all programs are chained
        while !ids.is_empty() {
            ids.retain(|&id| {
                let previous_id = programs[id].previous_id().unwrap();
                if let Some(mut chain) = chains.remove(previous_id) {
                    chain += 1;
                    chains.insert(id, chain);
                    false
                } else {
                    true
                }
            })
        }

        // return the appenadable PG ids
        match chains.len() {
            0 => Vec::new(),
            1 => vec![chains.iter().next().unwrap().0.to_string()],
            _ => {
                // if a len>1 chains exists return those, else return all
                if chains.values().any(|&l| l > 1) {
                    chains
                        .iter()
                        .filter_map(|(s, &l)| if l > 1 { Some(s.to_string()) } else { None })
                        .collect()
                } else {
                    chains.keys().map(|s| s.to_string()).collect()
                }
            }
        }
    }

    fn next_rumidup_id(&self) -> String {
        let mut count = 0;
        let id = |ct| {
            if ct > 0 {
                format!("rumidup.{ct}")
            } else {
                String::from("rumidup")
            }
        };

        let mut name = id(count);
        while self.0.programs().contains_key(&name) {
            name = id(count);
            count += 1;
        }

        name
    }
}

trait ProgramExt {
    fn is_markdup(&self) -> bool;
}

impl ProgramExt for Map<Program> {
    /// Currently detected software
    ///
    /// Picard MarkDups and variants
    /// sambamba markdup
    /// samtools markdup
    /// rumidup
    fn is_markdup(&self) -> bool {
        if let Some(name) = self.name() {
            if name.contains("MarkDuplicates") || name == "rumidup" {
                return true;
            }
            if let Some(command_line) = self.command_line() {
                if command_line.starts_with("samtools markdup") && name == "samtools" {
                    return true;
                }
            }
        } else if let Some(command_line) = self.command_line() {
            // sambamba doesn't set program name
            if command_line.starts_with("markdup") {
                return true;
            }
        }

        false
    }
}

#[cfg(test)]
mod test {
    use super::*;

    //use noodles_sam::header::{Header, record::value::map::{Map, Program}};

    const HEADER: &str = "@HD\tVN:1.0\tSO:coordinate
@SQ\tSN:1\tLN:248956422\tM5:2648ae1bacce4ec4b6cf337dcae37816\tUR:/tmp/Homo_sapiens.GRCh38.102/Homo_sapiens.GRCh38.dna.primary_assembly.fa
@RG\tID:H7L25DRX3.1\tCN:NKICMF\tLB:test_library\tPL:ILLUMINA\tSM:test_sample
@PG\tID:hisat2\tPN:hisat2\tCL:\"/tmp/envs/hisat2/bin/hisat2-align-s --wrapper basic-0 --rg CN:NKICMF --rg PL:ILLUMINA -x /tmp/grch38_snp_tran/genome_snp_tran --min-intronlen 20 --max-intronlen 500000 --rna-strandness FR -k 5 --minins 0 --rg-id H7L25DRX3.1 --rg LB:test_library --rg SM:test_sample --maxins 500 --fr --new-summary --threads 16 -1 /tmp/2130886.inpipe1 -2 /tmp/2130886.inpipe2\"\tVN:2.1.0
@PG\tID:samtools\tPN:samtools\tCL:samtools fixmate -u -m - -\tPP:hisat2\tVN:1.17
@PG\tID:samtools.1\tPN:samtools\tCL:samtools sort -m 3G -u -T /tmp/tmp.sL4pL151Jv -\tPP:samtools\tVN:1.17
@PG\tID:samtools.2\tPN:samtools\tCL:samtools view -@ 2 -C -T /tmp/Homo_sapiens.GRCh38.102/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o test_file.cram -\tPP:samtools.1\tVN:1.17
@PG\tID:samtools.3\tPN:samtools\tCL:samtools view -N /tmp/actreads.txt -o /tmp/actb_rna.bam test_file.cram\tPP:samtools.2\tVN:1.21
";

    const PG_RUMIDUP: &str = "@PG\tID:rumidup\tPN:rumidup\tCL:target/release/rumidup -b /tmp/actb_rna.bam -o /tmp/actb_rna_mdup.bam -p 2500\tPP:samtools.3\tVN:0.2.2";
    const PG_SAMBAMBA: &str = "@PG\tID:sambamba\tCL:markdup /tmp/actb_rna.bam /tmp/actb_rna_sambamba.bam\tPP:samtools.3\tVN:1.0";
    const PG_PICARD: &str = "@PG\tID:MarkDuplicates\tVN:Version:4.6.1.0\tCL:MarkDuplicates --INPUT /tmp/actb_rna.bam --OUTPUT /tmp/actb_rna_gatk.bam --METRICS_FILE /tmp/met.txt\tPN:MarkDuplicates";

    #[test]
    fn parse_and_add() {
        let mut h: BamHeader = HEADER.parse().unwrap();

        assert_eq!(h.as_ref().programs().iter().count(), 5);

        h.add_rumidup_pg("rumidup", "1.0");
        assert_eq!(h.as_ref().programs().iter().count(), 6);
        assert_eq!(h.as_ref().programs().last().unwrap().1.previous_id(), Some("samtools.3"));

        h.add_rumidup_pg("rumidup", "1.0");
        assert_eq!(h.as_ref().programs().last().unwrap().0, "rumidup.1");
        assert_eq!(h.as_ref().programs().last().unwrap().1.previous_id(), Some("rumidup"));
    }


    #[test]
    fn detect_previous_mdup() {
        let h: BamHeader = HEADER.parse().unwrap();
        assert!(h.detect_markdups().is_none());

        let pd = format!("{}{}", HEADER, PG_RUMIDUP);
        let h: BamHeader = pd.parse().unwrap();
        assert!(h.detect_markdups().is_some());

        let pd = format!("{}{}", HEADER, PG_SAMBAMBA);
        let h: BamHeader = pd.parse().unwrap();
        assert!(h.detect_markdups().is_some());

        let pd = format!("{}{}", HEADER, PG_PICARD);
        let h: BamHeader = pd.parse().unwrap();
        assert!(h.detect_markdups().is_some());
    }

    #[test]
    fn pg_chains() {
        let header = BamHeader::from(
            noodles_sam::Header::builder()
                .set_header(Default::default())
                .add_program("bwa", Map::<Program>::default())
                .add_program(
                    "samtools",
                    Map::<Program>::builder()
                    .set_previous_id("bwa")
                    .build()
                    .unwrap(),
                )
                .add_program(
                    "second_child",
                    Map::<Program>::builder()
                    .set_previous_id("second_root")
                    .build()
                    .unwrap(),
                )
                .add_program("second_root", Map::<Program>::default())
                .add_program(
                    "wrong_parent",
                    Map::<Program>::builder()
                    .set_previous_id("notexist")
                    .build()
                    .unwrap(),
                )
                .build(),
        );

        let chains = header.pg_chains();
        assert!(chains.len() == 2);
        assert!(chains.contains(&"samtools".to_owned()));
        assert!(chains.contains(&"second_child".to_owned()));
    }
}
