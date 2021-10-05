struct DupTracker {
    known_dups: HashSet<Vec<u8>>,
}

impl DupTracker {
    fn new() -> DupTracker {
        DupTracker { known_dups: HashSet::new }
    }

    fn is_known_dup(r: &UmiRecord) -> bool {
        self.known_dups.contains(r.seq_name);
    }
}



struct 
