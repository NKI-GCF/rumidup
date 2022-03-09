use std::fmt::{self, Debug};

use ahash::AHashMap;

use crate::bktree::BkTree;
pub use crate::bktree::Dist;

use crate::record::FragmentCoord;

#[derive(Clone, Debug)]
pub enum MarkResult {
    Unknown(CoordPair),
    Unique(CorrectedUmi),
    Duplicate(CorrectedUmi),
    OpticalDuplicate(CorrectedUmi),
    Failed,
    Skipped,
    Unusable,
}

#[derive(Clone, Debug, Default)]
pub struct CorrectedUmi(Vec<u8>);

impl fmt::Display for CorrectedUmi {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "{}", std::string::String::from_utf8_lossy(&self.0))
    }
}

impl From<&[u8]> for CorrectedUmi {
    fn from(i: &[u8]) -> CorrectedUmi {
        CorrectedUmi(i.to_vec())
    }
}

impl CorrectedUmi {
    pub fn is_corrected(&self) -> bool {
        !self.0.is_empty()
    }

    pub fn umi(&self) -> Option<&[u8]> {
        if self.0.is_empty() {
            None
        } else {
            Some(&self.0)
        }
    }

    pub fn from_cmp(my_umi: &[u8], correct_umi: &[u8]) -> Self {
        if my_umi == correct_umi {
            CorrectedUmi(Vec::new())
        } else {
            CorrectedUmi(correct_umi.to_vec())
        }
    }
}

pub type CoordPair = (FragmentCoord, FragmentCoord);
#[derive(Default)]
pub struct UmiClusters<T: Dist + Debug>(AHashMap<CoordPair, BkTree<T>>);

impl<T> FromIterator<(CoordPair, T)> for UmiClusters<T>
where
    T: Debug + Dist,
{
    fn from_iter<I: IntoIterator<Item = (CoordPair, T)>>(iter: I) -> Self {
        let mut trees = AHashMap::new();

        for (coord_pair, umiindex) in iter {
            let tree = trees.entry(coord_pair).or_insert_with(BkTree::new);
            tree.insert(umiindex);
        }

        UmiClusters(trees)
    }
}

impl<T> UmiClusters<T>
where
    T: Dist + Debug,
{
    pub fn insert(&mut self, coord_pair: CoordPair, value: T) {
        let tree = self.0.entry(coord_pair).or_insert_with(BkTree::new);
        tree.insert(value);
    }

    pub fn mark_duplicates(
        &self,
        coord: CoordPair,
        max_distance: usize,
    ) -> impl Iterator<Item = Vec<&T>> + '_ {
        let tree = &self.0[&coord];
        tree.partitions(max_distance)
            .map(|clust| clust.into_iter().flat_map(|i| tree[i].iter()).collect())
    }

    pub fn coordinate_pairs(&self) -> impl Iterator<Item = CoordPair> + '_ {
        self.0.keys().copied()
    }
}
