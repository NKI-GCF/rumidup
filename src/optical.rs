#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Coord {
    x: i32,
    y: i32,
}

impl Coord {
    pub fn dist(&self, other: &Coord) -> f32 {
        (((self.x - other.x).pow(2) + (self.y - other.y).pow(2)) as f32).sqrt()
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Location {
    tile: Vec<u8>,
    coord: Coord,
}

impl Location {
    pub fn new(tile: Vec<u8>, x: i32, y: i32) -> Location {
        Location {
            tile,
            coord: Coord { x, y },
        }
    }
}

pub struct DuplicateClusters {
    clusters: Vec<(usize, Location)>,
    boundaries: Vec<usize>,
}

impl DuplicateClusters {
    pub fn new<I: IntoIterator<Item = Location>>(i: I) -> DuplicateClusters {
        let mut clusters: Vec<_> = i.into_iter().enumerate().collect();
        clusters.sort_by(|a, b| a.1.cmp(&b.1));
        let boundaries = match clusters.len() {
            0 => vec![0],
            1 => vec![0, 1],
            n => {
                let mut v = vec![0];
                let mut current_tile = &clusters[0].1.tile;
                for (i, loc) in clusters.iter().enumerate() {
                    if current_tile != &loc.1.tile {
                        v.push(i);
                        current_tile = &loc.1.tile;
                    }
                }
                v.push(n);
                v
            }
        };

        DuplicateClusters {
            clusters,
            boundaries,
        }
    }

    pub fn count_optical_dups(&self, maxdist: i32) -> usize {
        let mut dups = 0;
        for (&start, &end) in self.boundaries.iter().zip(self.boundaries.iter().skip(1)) {
            let coords = self.clusters[start..end].iter().map(|l| l.1.coord);
            let mut opticals = OpticalClusters::new(coords);
            opticals.cluster(maxdist);
            dups += opticals.count_optical_dups();
        }
        dups
    }

    pub fn optical_cluster_ids(&self, maxdist: i32) -> Vec<usize> {
        let mut result = vec![0; self.clusters.len()];
        let mut id_offset = 0;
        for (&start, &end) in self.boundaries.iter().zip(self.boundaries.iter().skip(1)) {
            let tile_clusters = &self.clusters[start..end];

            let mut opticals = OpticalClusters::new(tile_clusters.iter().map(|c| c.1.coord));
            opticals.cluster(maxdist);

            for (index, id) in opticals.optical_cluster_ids().into_iter().enumerate() {
                let ori_pos = tile_clusters[index].0;
                result[ori_pos] = id + id_offset;
            }

            id_offset = result[tile_clusters.last().unwrap().0] + 1;
        }

        result
    }
}

pub struct OpticalClusters {
    clusters: Vec<Coord>,
    neighbours: Vec<Vec<usize>>,
    in_cluster: Vec<usize>,
    clustered: bool,
}

impl OpticalClusters {
    pub fn new<I: IntoIterator<Item = Coord>>(i: I) -> OpticalClusters {
        let clusters: Vec<_> = i.into_iter().collect();
        let in_cluster: Vec<_> = (0..clusters.len()).collect();
        let neighbours = in_cluster.iter().map(|&i| vec![i]).collect();

        OpticalClusters {
            clusters,
            neighbours,
            in_cluster,
            clustered: false,
        }
    }

    pub fn cluster(&mut self, maxdist: i32) {
        let maxdist = maxdist as f32;
        for i1 in 0..(self.neighbours.len() - 1) {
            for i2 in (i1 + 1)..(self.neighbours.len()) {
                let cluster_1 = self.in_cluster[i1];
                let cluster_2 = self.in_cluster[i2];

                if cluster_1 == cluster_2 {
                    continue;
                }

                if self.clusters[i1].dist(&self.clusters[i2]) < maxdist {
                    let adj_c2 = std::mem::take(&mut self.neighbours[cluster_2]);
                    // update cluster id for cluster 2 members
                    // add to cluster 1
                    for pos in adj_c2 {
                        self.in_cluster[pos] = cluster_1;
                        self.neighbours[cluster_1].push(pos);
                    }
                }
            }
        }
        self.clustered = true;
    }

    pub fn count_optical_dups(&self) -> usize {
        self.neighbours
            .iter()
            .map(|l| match l.len() {
                0 => 0,
                l => l - 1,
            })
            .sum()
    }

    pub fn optical_cluster_ids(&self) -> Vec<usize> {
        let mut result = vec![0; self.clusters.len()];
        for (id, neighbours) in self.neighbours.iter().filter(|l| !l.is_empty()).enumerate() {
            for &i in neighbours {
                result[i] = id;
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optical_dups() {
        let list = vec![
            Coord { x: 100, y: 100 },
            Coord { x: 100, y: 101 },
            Coord { x: 200, y: 200 },
        ];
        let mut o = OpticalClusters::new(list.clone());
        o.cluster(2);
        assert_eq!(o.count_optical_dups(), 1);
        assert_eq!(o.optical_cluster_ids(), vec![0, 0, 1]);

        let mut o = OpticalClusters::new(list.clone());
        o.cluster(200);
        assert_eq!(o.count_optical_dups(), 2);
        assert_eq!(o.optical_cluster_ids(), vec![0, 0, 0]);

        let mut o = OpticalClusters::new(list.clone());
        o.cluster(0);
        assert_eq!(o.count_optical_dups(), 0);
        assert_eq!(o.optical_cluster_ids(), vec![0, 1, 2]);
    }

    #[test]
    fn tile_duplicates() {
        let list = vec![
            Coord { x: 100, y: 100 },
            Coord { x: 100, y: 101 },
            Coord { x: 200, y: 200 },
        ];

        let onetile = vec![vec![1]; 3];

        let dc = DuplicateClusters::new(
            list.clone()
                .into_iter()
                .zip(onetile.into_iter())
                .map(|(coord, tile)| Location { tile, coord }),
        );

        assert_eq!(dc.count_optical_dups(50), 1);
        assert_eq!(dc.optical_cluster_ids(50), vec![0, 0, 1]);

        let twotiles = vec![vec![1], vec![2], vec![1]];
        let dc = DuplicateClusters::new(
            list.clone()
                .into_iter()
                .zip(twotiles.into_iter())
                .map(|(coord, tile)| Location { tile, coord }),
        );

        assert_eq!(dc.count_optical_dups(50), 0);
        assert_eq!(dc.count_optical_dups(150), 1);
        assert_eq!(dc.optical_cluster_ids(150), vec![0, 1, 0]);
    }
}
