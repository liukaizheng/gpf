use hashbrown::HashMap;

pub struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<usize>,
    pub n_groups: usize,
}

impl DisjointSet {
    pub fn new(n: usize) -> Self {
        Self {
            parent: Vec::from_iter(0..n),
            rank: vec![0; n],
            n_groups: n,
        }
    }

    fn link(&mut self, x: usize, y: usize) {
        if x == y {
            return;
        }

        if self.rank[x] > self.rank[y] {
            self.parent[y] = x;
        } else {
            self.parent[x] = y;
            if self.rank[x] == self.rank[y] {
                self.rank[y] += 1;
            }
        }
        self.n_groups -= 1;
    }

    #[inline(always)]
    pub fn merge(&mut self, x: usize, y: usize) {
        let x = self.find_set(x);
        let y = self.find_set(y);
        self.link(x, y);
    }

    fn find_set(&mut self, x: usize) -> usize {
        if x != self.parent[x] {
            self.parent[x] = self.find_set(self.parent[x]);
        }
        return self.parent[x];
    }

    pub fn output(&mut self) -> HashMap<usize, Vec<usize>> {
        let mut result = HashMap::<usize, Vec<usize>>::new();
        for i in 0..self.parent.len() {
            let p = self.find_set(i);
            result.entry(p).or_insert(vec![]).push(i);
        }
        result
    }
}
