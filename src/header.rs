use std::str::FromStr;

use ahash::AHashMap;
pub use noodles_sam::header::ParseError;
use noodles_sam::header::{Header, Program};

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

    pub fn add_rumidup_pg(&mut self) {
        let pps = self.pg_chains();

        if pps.is_empty() {
            let id = self.next_rumidup_id();
            self.0.programs_mut().insert(
                id.clone(),
                Program::builder()
                    .set_id(id)
                    .set_name("rumidup")
                    .set_version("0.1")
                    .set_command_line("cmd")
                    .build()
                    .expect("Error building program"),
            );
        } else {
            for p in pps {
                let id = self.next_rumidup_id();
                self.0.programs_mut().insert(
                    id.clone(),
                    Program::builder()
                        .set_id(id)
                        .set_name("rumidup")
                        .set_version("0.1")
                        .set_command_line("cmd")
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

impl ProgramExt for Program {
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
            if self.id().starts_with("sambamba") && command_line.starts_with("markdup") {
                return true;
            }
        }

        false
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use noodles_sam::header::{Header, Program};

    #[test]
    fn pg_chains() {
        let header = BamHeader::from(
            Header::builder()
                .set_header(noodles_sam::header::header::Header::default())
                .add_program(Program::builder().set_id("bwa").build().unwrap())
                .add_program(
                    Program::builder()
                        .set_id("samtools")
                        .set_previous_id("bwa")
                        .build()
                        .unwrap(),
                )
                .add_program(
                    Program::builder()
                        .set_id("second_child")
                        .set_previous_id("second_root")
                        .build()
                        .unwrap(),
                )
                .add_program(Program::builder().set_id("second_root").build().unwrap())
                .add_program(
                    Program::builder()
                        .set_id("wrong_parent")
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
