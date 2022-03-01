use std::convert::TryFrom;

use noodles_bam::record::{
    cigar::Cigar,
    data::{field::Value, Data, Field},
    Record as BamRecord,
};
use noodles_sam::record::{data::field::Tag, Cigar as SamCigar, Flags};
pub use noodles_sam::record::Position;
use thiserror::Error;

use crate::bktree::Dist;
use crate::optical::Location;

#[derive(Debug)]
pub struct UmiRecord {
    pub record: BamRecord,
    pub umi: Vec<u8>,
    pub mate_score: Option<u32>,
    pub mate_cigar: Option<Cigar>,
    pub location: Option<Location>,
}

#[derive(Debug, Default)]
struct FieldIndex {
    umi: Option<usize>,
    mate_score: Option<usize>,
    mate_cigar: Option<usize>,
}

impl TryFrom<&[u8]> for Location {
    type Error = UmiError;
    fn try_from(r: &[u8]) -> Result<Self, Self::Error> {
        // A01260:10:HWNYWDRXX:1:1273:8205:25598
        let mut e = r.split(|&b| b == b':')
            .skip(3)
            .take(4);

        let mut lanetile: Vec<u8> = e.next()
            .ok_or(UmiError::NoCoords)?
            .to_vec();
        lanetile.extend(e.next()
            .ok_or(UmiError::NoCoords)?
            .iter()
            .copied()
        );

        let x = String::from_utf8_lossy(e.next()
            .ok_or(UmiError::NoCoords)?)
            .parse().map_err(|_| UmiError::NoCoords)?;

        let y = String::from_utf8_lossy(e.next()
            .ok_or(UmiError::NoCoords)?)
            .parse().map_err(|_| UmiError::NoCoords)?;

        Ok(Location::new(lanetile, x, y))
    }
}

impl TryFrom<&Data> for FieldIndex {
    type Error = UmiError;
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

    pub fn read_name(&self) -> &[u8] {
        self.record.read_name()
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
            let len = self.record.cigar().reference_len().unwrap() as i32;
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
                let len = self.record.cigar().reference_len().unwrap() as i32;
                return Some(FragmentCoord::FragmentEnd(start + len));
            }
        }

        // two reads are mapped
        let mate_start = self.record.mate_position().map(|i| i.into())?;

        if flags.is_mate_reverse_complemented() {
            //maybe use template len for proper pairs
            let len = self.mate_cigar.as_ref()?.reference_len().ok()? as i32;
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

    pub fn score(&self) -> u32 {
        //FIXME cleanup?
        let rs = self
            .record
            .quality_scores()
            .scores()
            .map(|s| s.map(u8::from).unwrap_or(0) as u32)
            .sum();
        if let Some(ms) = self.mate_score {
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
        /*
        let data = self.record.data_mut();

        let old = data
            .insert(Field::new(Tag::UmiSequence, Value::String(tag.into())))
            .unwrap()
            .unwrap();
        data.insert(Field::new(
            Tag::OriginalUmiBarcodeSequence,
            Value::String(old.value().as_str().unwrap().to_owned()),
        ));
        */
    }
}

impl TryFrom<BamRecord> for UmiRecord {
    type Error = UmiError;

    fn try_from(mut r: BamRecord) -> Result<UmiRecord, Self::Error> {
        //extract data indices for required fields
        let fields = FieldIndex::try_from(r.data())?;
        let read_name = r.read_name();
        let location = Some(Location::try_from(read_name)?);

        let umi = if let Some(i) = fields.umi {
            r.data()
                .get_index(i)
                .unwrap()?
                .value()
                .as_str()
                .unwrap()
                .as_bytes()
                .to_vec()
        } else if let Some(umi) = umi_from_readname(read_name) {
            /*
            r.data_mut().insert(Field::new(
                Tag::UmiSequence,
                Value::String(String::from_utf8_lossy(&umi).to_string()),
            ));
            */
            umi
        } else {
            return Err(UmiError::NoUmi);
        };

        let mate_cigar = fields
            .mate_cigar
            .and_then(|i| r.data().get_index(i))
            .transpose()?
            .map(|value| {
                value
                    .value()
                    .as_str()
                    .ok_or(UmiError::InvalidMateCigar)
                    .and_then(|cs| {
                        cs.parse::<SamCigar>()
                            .map_err(|_| UmiError::InvalidMateCigar)
                    })
            })
            .transpose()?
            .map(|cigar| {
                let mut c = Vec::with_capacity(cigar.len());
                for &op in cigar.iter() {
                    let len = op.len();
                    let kind = op.kind() as u32;
                    let value = len << 4 | kind;
                    c.push(value);
                }
                Cigar::from(c)
            });

        let mate_score = fields
            .mate_score
            .and_then(|i| r.data().get_index(i))
            .transpose()?
            .and_then(|value| value.value().as_int())
            .map(|score| u32::try_from(score).map_err(|_| UmiError::NoMateScore))
            .transpose()?;


        Ok(UmiRecord {
            record: r,
            umi,
            mate_score,
            mate_cigar,
            location,
        })
    }
}

pub fn isbase(b: &u8) -> bool {
    b == &b'A' || b == &b'C' || b == &b'G' || b == &b'T' || b == &b'N'
}

fn umi_from_readname(r: &[u8]) -> Option<Vec<u8>> {
    r.split(|&c| c == b':').next_back().and_then(|maybe| {
        if maybe.len() > 3 && maybe.iter().all(isbase) {
            Some(maybe.to_vec())
        } else {
            None
        }
    })
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
pub enum UmiError {
    #[error("Error reading BAM")]
    IoError(#[from] std::io::Error),
    #[error("No MS tag in record data. Run samtools fixmate")]
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
}
