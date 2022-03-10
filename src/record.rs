use std::convert::TryFrom;

use noodles_bam::record::data::{field::Value, Data, Field};
pub use noodles_bam::record::Record as BamRecord;
use noodles_sam::record::{data::field::Tag, Cigar as SamCigar};
pub use noodles_sam::record::{Flags, Position};
use thiserror::Error;

use crate::bktree::Dist;
use crate::optical::Location;

pub type ReadName = Vec<u8>;

#[derive(Debug)]
pub struct UmiRecord {
    pub record: BamRecord,
    pub umi: Vec<u8>,
    pub mate_score: Option<u32>,
    pub mate_cigar: Option<SamCigar>,
    pub location: Option<Location>,
}

#[derive(Debug, Default)]
struct FieldIndex {
    umi: Option<usize>,
    mate_score: Option<usize>,
    mate_cigar: Option<usize>,
}

impl TryFrom<&Data> for FieldIndex {
    type Error = RecordError;
    fn try_from(data: &Data) -> Result<Self, Self::Error> {
        let mut fields = FieldIndex::default();
        for (index, key) in data.keys().enumerate() {
            let key = key?;
            match key {
                Tag::UmiSequence => fields.umi = Some(index),
                Tag::Other(_) if key.as_ref() == b"ms" => fields.mate_score = Some(index),
                Tag::MateCigar => fields.mate_cigar = Some(index),
                _ => {}
            }
        }
        Ok(fields)
    }
}

impl From<UmiRecord> for BamRecord {
    fn from(r: UmiRecord) -> BamRecord {
        r.record
    }
}

impl TryFrom<&[u8]> for Location {
    type Error = RecordError;
    fn try_from(name: &[u8]) -> Result<Location, RecordError> {
        // A01260:10:HWNYWDRXX:1:1273:8205:25598
        let mut e = name.split(|&b| b == b':').skip(3).take(4);

        let mut lanetile: Vec<u8> = e.next().ok_or(RecordError::NoCoords)?.to_vec();
        lanetile.extend(e.next().ok_or(RecordError::NoCoords)?.iter().copied());

        let x = String::from_utf8_lossy(e.next().ok_or(RecordError::NoCoords)?)
            .parse()
            .map_err(|_| RecordError::NoCoords)?;

        let y = String::from_utf8_lossy(e.next().ok_or(RecordError::NoCoords)?)
            .parse()
            .map_err(|_| RecordError::NoCoords)?;

        Ok(Location::new(lanetile, x, y))
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum FragmentCoord {
    Read1Start(i32),
    Read2Start(i32),
    MateStartFw(i32),
    MateStartRev(i32),
    FragmentEnd(i32),
}

impl FieldIndex {
    fn mate_score(&self, data: &Data) -> Result<Option<u32>, RecordError> {
        self.mate_score
            .and_then(|i| data.get_index(i))
            .transpose()?
            .and_then(|value| value.value().as_int())
            .map(|score| u32::try_from(score).map_err(|_| RecordError::NoMateScore))
            .transpose()
    }

    fn mate_cigar(&self, data: &Data) -> Result<Option<SamCigar>, RecordError> {
        self.mate_cigar
            .and_then(|i| data.get_index(i))
            .transpose()?
            .map(|value| {
                value
                    .value()
                    .as_str()
                    .ok_or(RecordError::InvalidMateCigar)
                    .and_then(|cs| {
                        cs.parse::<SamCigar>()
                            .map_err(|_| RecordError::InvalidMateCigar)
                    })
            })
            .transpose()
    }

    fn umi(&self, data: &Data) -> Result<Vec<u8>, RecordError> {
        if let Some(i) = self.umi {
            Ok(data
                .get_index(i)
                .unwrap()?
                .value()
                .as_str()
                .unwrap()
                .as_bytes()
                .to_vec())
        } else {
            Err(RecordError::NoUmi)
        }
    }

    fn set_umi(&self, umi: &[u8], data: &mut Data) -> Result<(), RecordError> {
        if data
            .insert(Field::new(
                Tag::UmiSequence,
                Value::String(String::from_utf8_lossy(umi).to_string()),
            ))
            .transpose()?
            .is_some()
        {
            return Err(RecordError::ParseError(
                "Record has an UMI both in tags and in the read name".to_string(),
            ));
        }
        Ok(())
    }
}

impl UmiRecord {
    pub fn from_bam_record(
        mut r: BamRecord,
        edit_readname: bool,
        extract_tags: bool,
        extract_location: bool,
    ) -> Result<UmiRecord, RecordError> {
        let data_fields = FieldIndex::try_from(r.data())?;
        // Assume the tag is in the readname. when found overwrite any tag in RX
        let umi = if let Some(umi) = umi_from_readname(r.read_name_mut(), edit_readname) {
            data_fields.set_umi(&umi, r.data_mut())?;
            umi
        } else {
            data_fields.umi(r.data())?
        };

        if umi.is_empty() {
            return Err(RecordError::NoUmi);
        }

        let mate_score = if extract_tags {
            Some(
                data_fields
                    .mate_score(r.data())?
                    .ok_or(RecordError::NoMateScore)?,
            )
        } else {
            None
        };

        let mate_cigar = if extract_tags {
            Some(
                data_fields
                    .mate_cigar(r.data())?
                    .ok_or(RecordError::NoMateCigar)?,
            )
        } else {
            None
        };

        let location = if extract_location {
            if Location::try_from(r.read_name()).is_err() {
                eprintln!("{:?}", r);
            }
            Some(Location::try_from(r.read_name())?)
        } else {
            None
        };

        Ok(UmiRecord {
            record: r,
            umi,
            mate_score,
            mate_cigar,
            location,
        })
    }

    pub fn position(&self) -> Option<Position> {
        self.record.position()
    }

    pub fn flags(&self) -> Flags {
        self.record.flags()
    }

    pub fn read_name(&self) -> &[u8] {
        self.record.read_name()
    }

    pub fn fragment_markers(&self) -> Option<(FragmentCoord, FragmentCoord)> {
        self.dedup_primary_coord().zip(self.dedup_secondary_coord())
    }

    pub fn is_paired(&self) -> bool {
        self.flags().is_segmented()
    }

    // determine the coord of the start of the read of the sequencing
    pub fn dedup_primary_coord(&self) -> Option<FragmentCoord> {
        let flags = self.record.flags();
        if flags.is_unmapped() || flags.is_supplementary() || flags.is_secondary() {
            return None;
        }

        let paired = flags.is_segmented();
        let start = self.position().map(|i| i.into())?;

        if flags.is_reverse_complemented() {
            let len = self.record.cigar().reference_len().unwrap() as i32;
            if !paired || flags.is_first_segment() {
                Some(FragmentCoord::Read1Start(start + len))
            } else {
                Some(FragmentCoord::Read2Start(start + len))
            }
        } else if !paired || flags.is_first_segment() {
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
        if !flags.is_segmented() || flags.is_mate_unmapped() {
            let start = self.record.position().map(|i| i.into())?;
            if flags.is_reverse_complemented() {
                return Some(FragmentCoord::FragmentEnd(start));
            } else {
                let len = self.record.cigar().reference_len().unwrap() as i32;
                return Some(FragmentCoord::FragmentEnd(start + len));
            }
        }

        // two reads are mapped
        let mate_start = self.record.mate_position().map(|i| i.into())?;

        if flags.is_mate_reverse_complemented() {
            //maybe use template len for proper pairs
            let len = self.mate_cigar.as_ref()?.reference_len() as i32;
            if flags.is_properly_aligned() {
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

    pub fn score(&self) -> u32 {
        self.record
            .quality_scores()
            .scores()
            .filter_map(|s| s.map(|q| u32::from(u8::from(q))).ok())
            .sum::<u32>()
            + self.mate_score.unwrap_or(0)
    }

    pub fn umi(&self) -> &[u8] {
        &self.umi
    }

    pub fn flag_dup(&mut self) {
        self.record.flags_mut().set(Flags::DUPLICATE, true);
    }

    pub fn unflag_dup(&mut self) {
        self.record.flags_mut().set(Flags::DUPLICATE, false);
    }

    pub fn correct_barcode_umi<T: Into<String>>(&mut self, umi: T, original_tag: bool) {
        let data = self.record.data_mut();

        let old = data
            .insert(Field::new(Tag::UmiSequence, Value::String(umi.into())))
            .unwrap()
            .unwrap();
        if original_tag {
            data.insert(Field::new(
                Tag::OriginalUmiBarcodeSequence,
                Value::String(old.value().as_str().unwrap().to_owned()),
            ));
        }
    }
}

pub fn isbase(b: &u8) -> bool {
    b == &b'A' || b == &b'C' || b == &b'G' || b == &b'T' || b == &b'N'
}

fn umi_from_readname(r: &mut Vec<u8>, edit: bool) -> Option<Vec<u8>> {
    if let Some(lastcolon) = r.iter().rev().position(|&c| c == b':') {
        let pos = r.len() - lastcolon;
        if lastcolon >= 3 && r[pos..].iter().all(isbase) {
            let umi = r[pos..].to_vec();
            if edit {
                r.truncate(pos - 1);
            }
            return Some(umi);
        }
    }

    None
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

#[derive(Debug, Error)]
pub enum RecordError {
    #[error("Error reading BAM")]
    IoError(#[from] std::io::Error),
    #[error("No MC (mate cigar) tag in record data. Run samtools fixmate")]
    NoMateCigar,
    #[error("The MC tag cannot be parsed")]
    InvalidMateCigar,
    #[error("No ms (mate score) tag in record data. Run samtools fixmate")]
    NoMateScore,
    #[error("No umi in record")]
    NoUmi,
    #[error("No umi in record")]
    NoReadName(#[from] std::ffi::FromBytesWithNulError),
    #[error("Error parsing coords from readname")]
    NoCoords,
    #[error("ParseError reading BAM: {0}")]
    ParseError(String),
}
