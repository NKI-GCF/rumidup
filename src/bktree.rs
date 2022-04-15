use std::ops::Index;

use smallvec::{smallvec, SmallVec};
use vec_map::VecMap;

pub trait Dist<Rhs = Self> {
    fn dist(&self, other: &Rhs) -> usize;
}

impl Dist for &str {
    fn dist(&self, b: &Self) -> usize {
        let la = self.len();
        let lb = b.len();
        if la == lb {
            self.chars()
                .zip(b.chars())
                .filter(|(a, b)| a == &'N' || b == &'N' || a != b)
                .count()
        } else {
            std::cmp::max(la, lb)
        }
    }
}

#[derive(Debug)]
pub struct UmiNode<T> {
    values: SmallVec<[T; 4]>,
    count: usize,
    edges: VecMap<usize>,
}

impl<T> UmiNode<T> {
    pub fn new(umi: T) -> UmiNode<T> {
        UmiNode {
            values: smallvec![umi],
            count: 1,
            edges: VecMap::new(),
        }
    }

    pub fn count(&self) -> usize {
        self.count
    }

    pub fn values(&self) -> &SmallVec<[T; 4]> {
        &self.values
    }
}

#[derive(Debug, Default)]
pub struct BkTree<T> {
    pub umis: Vec<UmiNode<T>>,
}

impl<T> BkTree<T>
where
    T: Dist,
{
    pub fn new() -> BkTree<T> {
        BkTree { umis: Vec::new() }
    }

    /// Insert a new umi into the `BkTree`. All inserted values with distance 0 will be stored in
    /// a `SmallVec` in the same node in the tree.
    ///
    /// Returns: boolean `true` when the umi not yet in the tree, `false` otherwise
    pub fn insert(&mut self, umi: T) -> bool {
        let pos = self.umis.len();
        if self.umis.is_empty() {
            self.umis.push(UmiNode::new(umi));
        } else {
            let mut current = 0;
            loop {
                if let Some(node) = self.umis.get_mut(current) {
                    let d = node.values[0].dist(&umi);
                    if d == 0 {
                        //already in tree
                        node.count += 1;
                        node.values.push(umi);
                        return false;
                    }
                    if let Some(target) = node.edges.get_mut(d) {
                        current = *target;
                    } else {
                        node.edges.insert(d, pos);
                        break;
                    }
                }
            }
            self.umis.push(UmiNode::new(umi));
        }
        true
    }

    /// Partition the UMI values in the tree grouping UMI up to `max_dist`.
    ///
    /// The UMI sequences are sorted (stable) based on number of occurences. Then, in order, a UMI is picked
    /// and all related UMI sequences up to `max_dist` are clusterd.
    ///
    /// Returns a `TreePartitions` iterator with the partitions.
    ///
    /// # Iterator behaviour
    /// Each iteration will return a `Vec<usize>`. The usize can be used to index the tree and and
    /// access the stored values. The first entry of the vec contains the most occurring UMI
    /// sequence. Ties are return in insertion order for reproducability.
    pub fn partitions(&self, max_dist: usize) -> TreePartitions<'_, T> {
        let mut node_order: Vec<_> = (0..self.umis.len()).collect();
        node_order.sort_by_key(|&e| self.umis[e].count);

        let work_edges_out: Vec<_> = self.umis.iter().map(|u| u.edges.clone()).collect();
        //let mut work_parents = vec![None; work_edges_out.len()];
        //work_edges_out.iter().enumerate().for_each(|(i, edges)| {
        //    edges.iter().for_each(|(dist, &to)| work_parents[to] = Some((i, dist)));
        //});

        TreePartitions {
            bktree: self,
            max_dist,
            partition: 0,
            node_order,
            assigned_partition: vec![0; self.umis.len()],
            edge_queue: Vec::new(),
            work_edges_out,
            //work_parents
        }
    }

    /// Return a `Vec<usize>` to children of the given node
    pub fn get_children(&self, start: usize) -> Vec<usize> {
        let mut children = Vec::new();
        let mut queue: Vec<usize> = self.umis[start].edges.values().copied().collect();
        while let Some(next) = queue.pop() {
            children.push(next);
            queue.extend(self.umis[next].edges.values().copied());
        }
        children
    }

    /// drain the tree returning the inserted values
    pub fn drain(self) -> impl Iterator<Item = T> {
        self.umis
            .into_iter()
            .flat_map(|node| node.values.into_iter())
    }
}

impl<T> Index<usize> for BkTree<T> {
    type Output = SmallVec<[T; 4]>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.umis[index].values
    }
}

/// Iterator over the partitions in the `BkTree`, created with the method `BkTree::partitions`
pub struct TreePartitions<'a, T> {
    bktree: &'a BkTree<T>,
    max_dist: usize,
    partition: usize,
    node_order: Vec<usize>,
    assigned_partition: Vec<usize>,
    edge_queue: Vec<usize>,
    work_edges_out: Vec<VecMap<usize>>,
    //work_parents: Vec<Option<(usize, usize)>>,
}

impl<'a, T> Iterator for TreePartitions<'a, T>
where
    T: Dist,
{
    type Item = Vec<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        let next_node = self.node_order.pop()?;
        if self.assigned_partition[next_node] != 0 {
            return self.next();
        }

        self.partition += 1;
        let max_dist = self.max_dist;
        let mut partition = Vec::new();

        let node = &self.bktree.umis[next_node];
        self.assigned_partition[next_node] = self.partition;
        partition.push(next_node);

        self.edge_queue.clear();
        self.edge_queue.push(0);

        while let Some(current) = self.edge_queue.pop() {
            let current_node = &self.bktree.umis[current];
            let d = node.values[0].dist(&current_node.values[0]);

            if d <= self.max_dist && self.assigned_partition[current] == 0 {
                //add this node to the partition
                self.assigned_partition[current] = self.partition;
                partition.push(current);
            }
            // follow the edges that might lead to other hits
            self.edge_queue.extend(
                self.work_edges_out[current]
                    .iter()
                    .filter(|&(dist, _)| abs_diff(dist, d) <= max_dist)
                    .map(|(_, target)| target),
            );
        }

        if partition.is_empty() {
            None
        } else {
            Some(partition)
        }
    }
}

impl<'a, T> TreePartitions<'a, T> {
    // FIXME tree pruning was a thing during development, code copied here but not fixed for
    // iterator, when performance is poor on large trees this may help
    #[allow(dead_code)]
    fn prune_tree(&mut self) {
        /*
        loop {
            if outgoing[bubble].is_empty() && bubble != 0 && result[bubble] != 0 {
                //eprintln!("leaf at {}", bubble);
                let (parent, dist) = parents[bubble].take().unwrap();
                outgoing[parent].remove(dist);
                bubble = parent;
            } else if outgoing[bubble].len() == 1 && bubble != 0 && result[bubble] != 0 {
                //eprintln!("single at {}", bubble);
                //remove my outgoing edge
                let (_dist, out) = outgoing[bubble].drain().next().unwrap();

                //get info from parent ege to bubble
                let (parent, dist) = parents[bubble].take().unwrap();

                //and modify on the parent
                outgoing[parent][dist] = out;
                parents[out] = Some((parent, dist));
                bubble = parent;
            } else {
                //eprintln!("bubble end");
                break;
            }

        }
        */
    }
}

#[inline]
fn abs_diff<T: Ord + std::ops::Sub<Output = T>>(a: T, b: T) -> T {
    if a < b {
        b - a
    } else {
        a - b
    }
}

pub fn dist(a: &[u8], b: &[u8]) -> usize {
    let la = a.len();
    let lb = b.len();
    if la == lb {
        a.iter()
            .zip(b.iter())
            .filter(|(&a, &b)| a == b'N' || b == b'N' || a != b)
            .count()
    } else {
        std::cmp::max(la, lb)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn it_works() {
        let mut tree = BkTree::new();
        tree.insert("CGAT");
        tree.insert("CGAT");
        tree.insert("CAAT");
        tree.insert("TTTA");
        println!("{:#?}", tree);
        assert_eq!(tree.umis.len(), 3);
    }

    #[test]
    fn children() {
        let mut tree = BkTree::new();
        tree.insert("CGCT");
        tree.insert("CCCT");
        tree.insert("CCAT");
        tree.insert("CGAT");

        let mut from_root = tree.get_children(0);
        from_root.sort();
        assert_eq!(from_root, vec![1, 2, 3]);
        assert_eq!(tree.get_children(2), vec![]);
        assert_eq!(tree.get_children(1), vec![3]);
    }

    #[test]
    fn partition2() {
        let mut tree = BkTree::new();
        tree.insert("CGCT"); // 0
        tree.insert("CCCT"); //1
        tree.insert("CCAT"); //2
        tree.insert("CGAT"); //3
        tree.insert("CGAT"); //-
        println!("{:#?}", tree);
        let mut d1 = tree.partitions(1);
        assert_eq!(d1.next(), Some(vec![3, 0, 2]));
        assert_eq!(d1.next(), Some(vec![1]));
        assert!(d1.next().is_none());

        let mut d2 = tree.partitions(2);
        assert_eq!(d2.next(), Some(vec![3, 0, 2, 1]));
        assert!(d2.next().is_none());
    }
}
