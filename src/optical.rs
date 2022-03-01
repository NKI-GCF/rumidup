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
        Location { tile, coord: Coord { x, y } }
    }
}

pub struct DuplicateClusters {
    clusters: Vec<Location>,
    boundaries: Vec<usize>,
}

impl DuplicateClusters {
    pub fn new<I: IntoIterator<Item=Location>>(i: I) -> DuplicateClusters {
        let mut clusters: Vec<_> = i.into_iter().collect();
        clusters.sort();
        let boundaries = match clusters.len() {
            0 => vec![0],
            1 => vec![0, 1],
            n => {
                let mut v = vec![0];
                let mut current_tile = &clusters[0].tile;
                for (i, loc) in clusters.iter().enumerate() {
                    if current_tile != &loc.tile {
                        v.push(i);
                        current_tile = &loc.tile;
                    }
                }
                v.push(n);
                v
            },
        };

        DuplicateClusters { clusters, boundaries }
    }

    pub fn count_optical_dups(&self, maxdist: i32) -> usize {
        let mut dups = 0;
        for (&start, &end) in self.boundaries.iter().zip(self.boundaries.iter().skip(1)) {
            let coords = self.clusters[start..end].iter().map(|l| l.coord);
            let mut opticals = OpticalClusters::new(coords);
            opticals.cluster(maxdist);
            dups += opticals.count_optical_dups();
        }
        dups
    }
}

pub struct OpticalClusters {
    clusters: Vec<Coord>,
    neighbours: Vec<Vec<usize>>,
    in_cluster: Vec<usize>,
    clustered: bool,
}

impl OpticalClusters {
    pub fn new<I: IntoIterator<Item=Coord>>(i: I) -> OpticalClusters {
        let clusters: Vec<_> = i.into_iter().collect();
        let in_cluster: Vec<_> = (0..clusters.len()).collect();
        let neighbours = in_cluster.iter().map(|&i| vec![i]).collect();

        OpticalClusters { clusters, neighbours, in_cluster, clustered: false }
    }

    pub fn cluster(&mut self, maxdist: i32) {
        let maxdist = maxdist as f32;
        for i1 in 0..(self.neighbours.len() - 1) {
            for i2 in (i1+1)..(self.neighbours.len()) {
                let cluster_1 = self.in_cluster[i1];
                let cluster_2 = self.in_cluster[i2];

                if cluster_1 == cluster_2 {
                    continue;
                }

                if self.clusters[i1].dist(&self.clusters[i2]) < maxdist {
                    let adj_c2 = std::mem::replace(&mut self.neighbours[cluster_2], Vec::new());
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
            .map(|l| { 
                match l.len() {
                    0 => 0,
                    l => l - 1
                }
            })
            .sum()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optical_dups() {
        let list = vec![Coord {x: 100, y: 100}, Coord {x: 100, y: 101}, Coord {x: 200, y: 200}];
        let mut o = OpticalClusters::new(list.clone());
        o.cluster(2);
        assert_eq!(o.count_optical_dups(), 1);

        let mut o = OpticalClusters::new(list.clone());
        o.cluster(200);
        assert_eq!(o.count_optical_dups(), 2);

        let mut o = OpticalClusters::new(list.clone());
        o.cluster(0);
        assert_eq!(o.count_optical_dups(), 0);
    }

}

