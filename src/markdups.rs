pub enum MarkResult {
    Unique,
    Duplicate(CorrectedUmi),
}

pub struct CorrectedUmi(Vec<u8>);


