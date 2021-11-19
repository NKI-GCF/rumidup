use std::convert::TryFrom;

use thiserror::Error;
use noodles_sam::record::{
    cigar::Cigar,
    data::field::tag::Tag,
    data::field::{Field, Value},
    position::Position,
    Flags, Record as SamRecord,
};
use umibk::Dist;

#[derive(Debug)]
pub struct UmiRecord {
    pub record: SamRecord,
    pub umi: Vec<u8>,
}

impl From<UmiRecord> for SamRecord {
    fn from(r: UmiRecord) -> SamRecord {
        r.record
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum FragmentCoord {
    Read1Start(i32),
    Read2Start(i32),
    MateStartFw(i32),
    MateStartRev(i32),
    FragmentEnd(i32),
}

impl UmiRecord {
    pub fn position(&self) -> Option<Position> {
        self.record.position()
    }

    pub fn flags(&self) -> Flags {
        self.record.flags()
    }

    pub fn read_name(&self) -> &str {
        self.record.read_name().map(|n| n.as_ref()).unwrap_or("")
    }

    pub fn fragment_markers(&self) -> Option<(FragmentCoord, FragmentCoord)> {
        self.dedup_primary_coord().zip(self.dedup_secondary_coord())
    }

    pub fn is_paired(&self) -> bool {
        self.flags().is_paired()
    }

    // determine the coord of the start of the read of the sequencing
    pub fn dedup_primary_coord(&self) -> Option<FragmentCoord> {
        let flags = self.record.flags();
        if flags.is_unmapped() || flags.is_supplementary() || flags.is_secondary() {
            return None;
        }

        let paired = flags.is_paired();
        let start = self.position().map(|i| i.into())?;

        if flags.is_reverse_complemented() {
            let len = self.record.cigar().reference_len() as i32;
            if !paired || flags.is_read_1() {
                Some(FragmentCoord::Read1Start(start + len))
            } else {
                Some(FragmentCoord::Read2Start(start + len))
            }
        } else if !paired || flags.is_read_1() {
            Some(FragmentCoord::Read1Start(start))
        } else {
            Some(FragmentCoord::Read2Start(start))
        }
    }

    //given a read figure out the mate coord that will be used to
    //deduplicate the read/readpair
    pub fn dedup_secondary_coord(&self) -> Option<FragmentCoord> {
        let flags = self.record.flags();

        //read is single end (no mate or unmapped)
        if !flags.is_paired() || flags.is_mate_unmapped() {
            let start = self.record.position().map(|i| i.into())?;
            if flags.is_reverse_complemented() {
                return Some(FragmentCoord::FragmentEnd(start));
            } else {
                let len = self.record.cigar().reference_len() as i32;
                return Some(FragmentCoord::FragmentEnd(start + len));
            }
        }

        // two reads are mapped
        let mate_start = self.record.mate_position().map(|i| i.into())?;

        if flags.is_mate_reverse_complemented() {
            //maybe use template len for proper pairs
            let len = self.mate_cigar()?.reference_len() as i32;
            if flags.is_proper_pair() {
                //test bam sanity
                let start: i32 = self.record.position().map(|i| i.into())?;
                debug_assert_eq!(mate_start + len - start, self.record.template_length())
            }
            Some(FragmentCoord::MateStartRev(mate_start + len))
        } else {
            Some(FragmentCoord::MateStartFw(mate_start))
        }
    }

    pub fn is_dedup_candidate(&self) -> bool {
        let flags = self.record.flags();
        !flags.is_unmapped() && !flags.is_supplementary() && !flags.is_secondary()
    }

    //FIXME use as_int
    pub fn ms(&self) -> Option<u32> {
        if self.record.flags().is_mate_unmapped() {
            None
        } else {
            self.record
                .data()
                .get(Tag::try_from(*b"ms").unwrap())
                .and_then(|f| f.value().as_int())
                .and_then(|v| u32::try_from(v).ok())
        }
    }

    pub fn mate_cigar(&self) -> Option<Cigar> {
        self.record
            .data()
            .get(Tag::MateCigar)
            .and_then(|cs| cs.value().as_str())
            .and_then(|s| s.parse().ok())
    }

    pub fn score(&self) -> u32 {
        //FIXME cleanup?
        let rs = self
            .record
            .quality_scores()
            .iter()
            .map(|&s| u8::from(s) as u32)
            .sum();
        if let Some(ms) = self.ms() {
            rs + ms
        } else {
            rs
        }
    }

    pub fn flag_dup(&mut self) {
        self.record.flags_mut().set(Flags::DUPLICATE, true);
    }

    pub fn unflag_dup(&mut self) {
        self.record.flags_mut().set(Flags::DUPLICATE, false);
    }

    pub fn correct_barcode_tag<T: Into<String>>(&mut self, tag: T) {
        let data = self.record.data_mut();

        let old = data.insert(Field::new(Tag::UmiSequence, Value::String(tag.into()))).unwrap();
        data.insert(
            Field::new(Tag::OriginalUmiBarcodeSequence, Value::String(old.value().as_str().unwrap().to_owned()))
        );
    }
}

impl TryFrom<SamRecord> for UmiRecord {
    type Error = UmiError;

    fn try_from(mut r: SamRecord) -> Result<UmiRecord, Self::Error> {
        // attempt 1 umi from read name
        let maybe_umi = r
            .read_name()
            .and_then(|name| name.split(':').next_back())
            .map(|v| v.as_bytes().to_vec());
        if let Some(maybe_umi) = maybe_umi {
            if maybe_umi.len() > 3 && maybe_umi.iter().all(isbase) {
                let umi = maybe_umi.to_vec();
                r.data_mut().insert(Field::new(Tag::UmiSequence, Value::String(String::from_utf8_lossy(&umi).to_string())));
                return Ok(UmiRecord {
                    record: r,
                    umi,
                });
            }
        }
        if let Some(umi) = record_umi_tag(&r) {
            return Ok(UmiRecord { record: r, umi });
        }

        Err(UmiError::NoUmi)
    }
}

impl Dist for UmiRecord {
    fn dist(&self, other: &UmiRecord) -> usize {
        self.umi
            .iter()
            .zip(other.umi.iter())
            .filter(|(a, b)| a != b)
            .count()
    }
}

pub fn isbase(b: &u8) -> bool {
    b == &b'A' || b == &b'C' || b == &b'G' || b == &b'T' || b == &b'N'
}

// create extension  trait for record?
pub fn record_umi_tag(record: &SamRecord) -> Option<Vec<u8>> {
    record
        .data()
        .get(Tag::UmiSequence)
        .and_then(|f| f.value().as_str())
        .map(|v| v.as_bytes().to_vec())
}

#[derive(Debug, Error)]
pub enum UmiError {
    #[error("No umi in record")]
    NoUmi,
    #[error("No umi in record")]
    NoReadName(#[from] std::ffi::FromBytesWithNulError),
}
