use noodles_sam::header::record::value::map::program::tag as program_tag;
pub use noodles_sam::header::ParseError;
use noodles_sam::header::{
    record::value::map::{Map, Program},
    Header,
};

pub struct BamHeader(Header);

impl From<Header> for BamHeader {
    fn from(h: Header) -> BamHeader {
        BamHeader(h)
    }
}

impl AsRef<Header> for BamHeader {
    fn as_ref(&self) -> &Header {
        &self.0
    }
}

impl BamHeader {
    // FIXME: Check debug formatting or write custom formatter
    pub fn detect_markdups(&self) -> Option<String> {
        self.0
            .programs()
            .as_ref()
            .values()
            .find(|p| p.is_markdup())
            .map(|p| format!("{:?}", p))
    }

    // FIXME convert to Result function
    pub fn add_rumidup_pg(&mut self, command_line: &str, version: &str) {
        self.0
            .programs_mut()
            .add(
                "rumidup",
                Map::<Program>::builder()
                    .insert(program_tag::NAME, b"rumidup")
                    .insert(program_tag::VERSION, version.as_bytes())
                    .insert(program_tag::COMMAND_LINE, command_line.as_bytes())
                    .build()
                    .expect("Error building program"),
            )
            .expect("Error adding program to header");
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
        if let Some(name) = self.other_fields().get(&program_tag::NAME) {
            if name == "MarkDuplicates" || name == b"rumidup" {
                return true;
            }
            if let Some(command_line) = self.other_fields().get(&program_tag::COMMAND_LINE) {
                if command_line.starts_with(b"samtools markdup") && name == b"samtools" {
                    return true;
                }
            }
        } else if let Some(command_line) = self.other_fields().get(&program_tag::COMMAND_LINE) {
            // sambamba doesn't set program name
            if command_line.starts_with(b"markdup") {
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
        let mut h: BamHeader = HEADER.parse::<Header>().unwrap().into();

        assert_eq!(h.as_ref().programs().as_ref().iter().count(), 5);

        h.add_rumidup_pg("rumidup", "1.0");
        assert_eq!(h.as_ref().programs().as_ref().iter().count(), 6);
        assert!(
            matches!(h.as_ref().programs().as_ref().last().unwrap().1.other_fields().get(&program_tag::PREVIOUS_PROGRAM_ID),
            Some(v) if v == b"samtools.3")
        );

        h.add_rumidup_pg("rumidup", "1.0");
        assert_eq!(
            h.as_ref().programs().as_ref().last().unwrap().0,
            "rumidup-rumidup"
        );
        assert!(
            matches!(h.as_ref().programs().as_ref().last().unwrap().1.other_fields().get(&program_tag::PREVIOUS_PROGRAM_ID),
            Some(v) if v == b"rumidup")
        );
    }

    #[test]
    fn detect_previous_mdup() {
        let h: BamHeader = HEADER.parse::<Header>().unwrap().into();
        assert!(h.detect_markdups().is_none());

        let pd = format!("{}{}", HEADER, PG_RUMIDUP);
        let h: BamHeader = pd.parse::<Header>().unwrap().into();
        assert!(h.detect_markdups().is_some());

        let pd = format!("{}{}", HEADER, PG_SAMBAMBA);
        let h: BamHeader = pd.parse::<Header>().unwrap().into();
        assert!(h.detect_markdups().is_some());

        let pd = format!("{}{}", HEADER, PG_PICARD);
        let h: BamHeader = pd.parse::<Header>().unwrap().into();
        assert!(h.detect_markdups().is_some());
    }
}
