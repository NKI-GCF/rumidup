pub use noodles_core::Position;
use noodles_sam::alignment::record::data::field::Tag;
pub use noodles_sam::alignment::record::Flags;
pub use noodles_sam::alignment::record_buf::RecordBuf as BamRecord;
use noodles_sam::alignment::record_buf::{data::field::Value, Cigar};
use noodles_sam::record::Cigar as RawCigar;
use thiserror::Error;

use crate::bktree::Dist;
use crate::optical::Location;

#[derive(Debug)]
pub struct UmiRecord {
    pub record: BamRecord,
    pub umi: Vec<u8>,
    pub mate_score: Option<i32>,
    pub mate_cigar: Option<Cigar>,
    pub location: Option<Location>,
}

impl From<UmiRecord> for BamRecord {
    fn from(r: UmiRecord) -> BamRecord {
        r.record
    }
}

/// Illumina readname to Location converter
impl TryFrom<&[u8]> for Location {
    type Error = RecordError;
    fn try_from(name: &[u8]) -> Result<Location, RecordError> {
        // A01260:10:HWNYWDRXX:1:1273:8205:25598
        let mut e = name.split(|&b| b == b':').skip(3).take(4);

        let mut lanetile: Vec<u8> = e.next().ok_or(RecordError::NoCoords)?.to_vec();
        lanetile.push(b':');
        lanetile.extend(e.next().ok_or(RecordError::NoCoords)?.iter().copied());

        let x: i32 = lexical_core::parse(e.next().ok_or(RecordError::NoCoords)?)
            .map_err(|_| RecordError::NoCoords)?;

        let y = lexical_core::parse(e.next().ok_or(RecordError::NoCoords)?)
            .map_err(|_| RecordError::NoCoords)?;

        Ok(Location::new(lanetile, x, y))
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum FragmentCoord {
    Read1Start(usize),
    Read2Start(usize),
    MateStartFw(usize),
    MateStartRev(usize),
    FragmentEnd(usize),
}

impl From<BamRecord> for UmiRecord {
    fn from(record: BamRecord) -> UmiRecord {
        UmiRecord {
            record,
            umi: Vec::with_capacity(10),
            mate_score: None,
            mate_cigar: None,
            location: None,
        }
    }
}

impl UmiRecord {
    /// Extract the UMI sequence from the record
    /// If the RX tag is present this is used. Otherwise the final part of the readname is
    /// considered to be the UMI sequence (this is the BCLconvert default). If this looks like a
    /// umi the sequence is put into an RX tag and (optionally) removed from the readname
    pub fn extract_umi(&mut self, edit_readname: bool) -> Result<(), RecordError> {
        if let Some(rx) = self.record.data().get(&Tag::UMI_SEQUENCE) {
            if let Value::String(ref umi) = rx {
                self.umi.extend(umi.as_slice());
            } else {
                return Err(RecordError::NoUmi);
            }
        } else if let Some(read_name) = self
            .record
            .name()
            .map(|r| AsRef::<[u8]>::as_ref(r).to_vec())
        {
            if let Some((umi, clipped)) = umi_from_readname(&read_name) {
                if edit_readname {
                    self.record.name_mut().replace(clipped.into());
                }
                self.umi.extend(umi);
                self.record
                    .data_mut()
                    .insert(Tag::UMI_SEQUENCE, Value::String(umi.to_vec().into()));
            } else {
                return Err(RecordError::NoUmi);
            }
        } else {
            return Err(RecordError::NoUmi);
        }

        Ok(())
    }

    /// Paired end reads require both MC and ms tags for rumidup to work
    pub fn extract_mate_tags(&mut self) -> Result<(), RecordError> {
        let data = self.record.data();

        let ms_value = data
            .get(&Tag::from(*b"ms"))
            .ok_or(RecordError::NoMateScore)?;
        let ms = ms_value
            .as_int()
            .ok_or(RecordError::MateScoreOutOfRange)
            .and_then(|ms| i32::try_from(ms).map_err(|_| RecordError::MateScoreOutOfRange))?;
        self.mate_score = Some(ms);

        self.mate_cigar = Some(
            data.get(&Tag::MATE_CIGAR)
                .ok_or(RecordError::NoMateCigar)
                .and_then(|v| {
                    if let Value::String(cigar) = v {
                        RawCigar::new(cigar.as_slice())
                            .try_into()
                            .map_err(|_| RecordError::NoMateCigar)
                    } else {
                        Err(RecordError::NoMateCigar)
                    }
                })?,
        );
        Ok(())
    }

    /// Extract the tile and coordinates when illumina read names are used. Can then be used for
    /// optical duplicate detection.
    pub fn extract_location(&mut self) -> Result<(), RecordError> {
        self.location = Some(Location::try_from(self.read_name())?);
        Ok(())
    }

    pub fn position(&self) -> Option<Position> {
        self.record.alignment_start()
    }

    pub fn flags(&self) -> Flags {
        self.record.flags()
    }

    pub fn read_name(&self) -> &[u8] {
        self.record.name().unwrap().as_ref()
    }

    pub fn display_name(&self) -> String {
        String::from_utf8_lossy(self.read_name()).to_string()
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
        let start: usize = self.position().map(|i| i.into())?;

        if flags.is_reverse_complemented() {
            let len = self.record.cigar().alignment_span();
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
            let start = self.record.alignment_start().map(|i| i.into())?;
            if flags.is_reverse_complemented() {
                return Some(FragmentCoord::FragmentEnd(start));
            } else {
                let len = self.record.cigar().alignment_span();
                return Some(FragmentCoord::FragmentEnd(start + len));
            }
        }

        // two reads are mapped
        let mate_start = self.record.mate_alignment_start().map(|i| i.into())?;

        if flags.is_mate_reverse_complemented() {
            //maybe use template len for proper pairs
            let len = self.mate_cigar.as_ref()?.alignment_span();
            Some(FragmentCoord::MateStartRev(mate_start + len))
        } else {
            Some(FragmentCoord::MateStartFw(mate_start))
        }
    }

    pub fn is_dedup_candidate(&self) -> bool {
        let flags = self.record.flags();
        !flags.is_unmapped() && !flags.is_supplementary() && !flags.is_secondary()
    }

    pub fn score(&self) -> i32 {
        self.record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&s| i32::from(s))
            .sum::<i32>()
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

    pub fn correct_barcode_umi(&mut self, umi: &str, original_tag: bool) {
        let data = self.record.data_mut();

        let old = data
            .insert(Tag::UMI_SEQUENCE, Value::String(umi.into()))
            .unwrap();
        if original_tag {
            data.insert(Tag::ORIGINAL_UMI_BARCODE_SEQUENCE, old.1);
        }
    }
}

pub fn isbase(b: &u8) -> bool {
    b == &b'A' || b == &b'C' || b == &b'G' || b == &b'T' || b == &b'N'
}

fn umi_from_readname(r: &[u8]) -> Option<(&[u8], &[u8])> {
    if let Some(lastcolon) = r.iter().rev().position(|&c| c == b':') {
        let pos = r.len() - lastcolon;
        if lastcolon >= 3
            && r[pos..]
                .iter()
                .filter(|&&c| c != b'+' && c != b'_')
                .all(isbase)
        {
            let umi = &r[pos..];
            let clipped = &r[0..pos - 1];
            return Some((umi, clipped));
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
    #[error("Tag ms (mate score) out of range")]
    MateScoreOutOfRange,
    #[error("No umi in record")]
    NoUmi,
    #[error("No readname in record")]
    NoReadName(#[from] std::ffi::FromBytesWithNulError),
    #[error("Error parsing coords from readname")]
    NoCoords,
    #[error("ParseError reading BAM: {0}")]
    ParseError(String),
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn umi_readname() {
        assert_eq!(
            umi_from_readname(&b"A01260:10:HWNYWDRXX:1:1273:8205:25598:CGACCTAGN"[..]),
            Some((
                &b"CGACCTAGN"[..],
                &b"A01260:10:HWNYWDRXX:1:1273:8205:25598"[..]
            ))
        );
        assert_eq!(
            umi_from_readname(&b"A01260:10:HWNYWDRXX:1:1273:8205:25198:CGACC+TAGNA"[..]),
            Some((
                &b"CGACC+TAGNA"[..],
                &b"A01260:10:HWNYWDRXX:1:1273:8205:25198"[..]
            ))
        );

        assert_eq!(
            umi_from_readname(&b"A01260:10:HWNYWDRXX:1:1273:8205:25198:CGACC_TAGNA"[..]),
            Some((
                &b"CGACC_TAGNA"[..],
                &b"A01260:10:HWNYWDRXX:1:1273:8205:25198"[..]
            ))
        );

        assert_eq!(
            umi_from_readname(&b"A01260:10:HWNYWDRXX:1:1273:8205:25598:CGACCXAGC"[..]),
            None
        );

        assert_eq!(
            umi_from_readname(&b"A01260:10:HWNYWDRXX:1:1273:8205:25598:CC"[..]),
            None
        );
    }

    #[test]
    fn location() {
        let loc =
            Location::try_from(&b"A01260:10:HWNYWDRXX:1:1273:8205:25198:CGACC+TAGNA"[..]).unwrap();
        let testloc = Location::new(b"1:1273".to_vec(), 8205, 25198);

        assert_eq!(loc, testloc);
    }
}
