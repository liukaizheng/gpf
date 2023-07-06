use crate::INVALID_IND;

#[derive(Clone)]
struct Arc {
    head: usize,
    next: usize,
    sister: usize,
    /// residual capacity
    r_cap: f64,
}

impl Default for Arc {
    fn default() -> Self {
        Self {
            head: INVALID_IND,
            next: INVALID_IND,
            sister: INVALID_IND,
            r_cap: 0.0,
        }
    }
}

impl Arc {
    fn new(head: usize, next: usize, sister: usize, r_cap: f64) -> Self {
        Self {
            head,
            next,
            sister,
            r_cap,
        }
    }
}

const ORPHAN: usize = INVALID_IND - 1;
const TERMINAL: usize = INVALID_IND - 2;

pub struct GraphCut {
    time: usize,
    queue_first: [usize; 2],
    queue_last: [usize; 2],
    orphan_first: usize,
    orphan_last: usize,
    pub is_sink: Vec<bool>,
    /// residual capacities of nodes
    tr_cap: Vec<f64>,
    flow: f64,
    /// all arcs
    arcs: Vec<Arc>,
    /// distance to the terminal
    dist: Vec<usize>,
    /// first arcs of nodes
    first_arc: Vec<usize>,
    /// next active node
    next: Vec<usize>,
    /// next orphan node
    next_orphan: Vec<usize>,
    /// parent arcs of nodes
    parent: Vec<usize>,
    /// time stamp showing when dist was computed
    ts: Vec<usize>,
}

impl GraphCut {
    pub fn new(source_cap: &[f64], sink_cap: &[f64]) -> Self {
        let n_nodes = source_cap.len();
        let dist = vec![1; n_nodes];
        let is_sink = vec![false; n_nodes];
        let next = vec![INVALID_IND; n_nodes];
        let next_orphan = vec![INVALID_IND; n_nodes];
        let parent = vec![TERMINAL; n_nodes];
        let ts = vec![0; n_nodes];
        let mut flow = 0.0;
        let tr_cap = Vec::from_iter(source_cap.iter().zip(sink_cap).map(|(&src, &sink)| {
            flow += if src < sink { src } else { sink };
            src - sink
        }));
        let mut graph = Self {
            time: 0,
            queue_first: [INVALID_IND; 2],
            queue_last: [INVALID_IND; 2],
            orphan_first: INVALID_IND,
            orphan_last: INVALID_IND,
            is_sink,
            tr_cap,
            flow,
            arcs: Vec::new(),
            dist,
            first_arc: vec![INVALID_IND; n_nodes],
            next,
            next_orphan,
            parent,
            ts,
        };

        for i in 0..n_nodes {
            graph.is_sink[i] = graph.tr_cap[i] < 0.0;
            if graph.tr_cap[i] == 0.0 {
                graph.parent[i] = INVALID_IND;
            } else {
                graph.set_active(i);
            }
        }
        graph
    }

    #[inline(always)]
    fn set_active(&mut self, i: usize) {
        if self.next[i] != INVALID_IND {
            return;
        }
        if self.queue_last[1] == INVALID_IND {
            self.queue_first[1] = i;
        } else {
            self.next[self.queue_last[1]] = i;
        }
        self.queue_last[1] = i;
        self.next[i] = i;
    }

    #[inline]
    fn next_active(&mut self) -> usize {
        loop {
            let mut i = self.queue_first[0];
            if i == INVALID_IND {
                i = self.queue_first[1];
                self.queue_first[0] = self.queue_first[1];
                self.queue_last[0] = self.queue_last[1];
                self.queue_first[1] = INVALID_IND;
                self.queue_last[1] = INVALID_IND;
                if i == INVALID_IND {
                    return INVALID_IND;
                }
            }
            if self.next[i] == i {
                self.queue_first[0] = INVALID_IND;
                self.queue_last[0] = INVALID_IND;
            } else {
                self.queue_first[0] = self.next[i];
            }
            self.next[i] = INVALID_IND;
            if self.parent[i] != INVALID_IND {
                return i;
            }
        }
    }

    #[inline(always)]
    fn set_orphan_front(&mut self, i: usize) {
        self.parent[i] = ORPHAN;
        self.next_orphan[i] = self.orphan_first;
        self.orphan_first = i;
    }

    #[inline(always)]
    fn set_orphan_rear(&mut self, i: usize) {
        self.parent[i] = ORPHAN;
        if self.orphan_last != INVALID_IND {
            self.next_orphan[self.orphan_last] = i;
        } else {
            self.orphan_first = i;
        }
        self.orphan_last = i;
        self.next_orphan[i] = INVALID_IND;
    }

    #[inline(always)]
    pub fn add_edge(&mut self, i: usize, j: usize, cap: f64, rev_cap: f64) {
        let arc_id = self.arcs.len();
        let sisiter_arc_id = arc_id + 1;
        self.arcs
            .push(Arc::new(j, self.first_arc[i], sisiter_arc_id, cap)); // arc
        self.first_arc[i] = arc_id;
        self.arcs
            .push(Arc::new(i, self.first_arc[j], arc_id, rev_cap)); // sisiter arc
        self.first_arc[j] = sisiter_arc_id;
    }

    fn process_source_orphan(&mut self, i: usize) {
        let mut a0_min = INVALID_IND;
        let mut d_min = INVALID_IND;
        let mut a0 = self.first_arc[i];
        while a0 != INVALID_IND {
            if self.arcs[self.arcs[a0].sister].r_cap == 0.0 {
                a0 = self.arcs[a0].next;
                continue;
            }
            let mut j = self.arcs[a0].head;
            if self.is_sink[j] || self.parent[j] == INVALID_IND {
                a0 = self.arcs[a0].next;
                continue;
            }
            // checking the origin of j
            let mut d = 0;
            loop {
                if self.ts[j] == self.time {
                    d += self.dist[j];
                    break;
                }
                let a = self.parent[j];
                d += 1;
                if a == TERMINAL {
                    self.ts[j] = self.time;
                    self.dist[j] = 1;
                    break;
                }
                if a == ORPHAN {
                    d = INVALID_IND;
                    break;
                }
                j = self.arcs[a].head;
            }
            if d != INVALID_IND {
                // j originates from the source - done
                if d < d_min {
                    a0_min = a0;
                    d_min = d;
                }
                // set marks along the path
                j = self.arcs[a0].head;
                loop {
                    if self.ts[j] == self.time {
                        break;
                    }
                    self.ts[j] = self.time;
                    self.dist[j] = d;
                    d -= 1;
                    j = self.arcs[self.parent[j]].head;
                }
            }
            a0 = self.arcs[a0].next;
        }

        self.parent[i] = a0_min;
        if self.parent[i] != INVALID_IND {
            self.ts[i] = self.time;
            self.dist[i] = d_min + 1;
        } else {
            // process neighbors
            let mut a0 = self.first_arc[i];
            while a0 != INVALID_IND {
                let j = self.arcs[a0].head;
                let a = self.parent[j];
                if !self.is_sink[j] && a != INVALID_IND {
                    if self.arcs[self.arcs[a0].sister].r_cap != 0.0 {
                        self.set_active(j);
                    }
                    if a != TERMINAL && a != ORPHAN && self.arcs[a].head == i {
                        self.set_orphan_rear(j); // add j to the end of the adoption list
                    }
                }
                a0 = self.arcs[a0].next;
            }
        }
    }

    fn process_sink_orphan(&mut self, i: usize) {
        let mut a0_min = INVALID_IND;
        let mut d_min = INVALID_IND;
        // trying to find a new parent
        let mut a0 = self.first_arc[i];
        while a0 != INVALID_IND {
            if self.arcs[a0].r_cap == 0.0 {
                a0 = self.arcs[a0].next;
                continue;
            }
            let mut j = self.arcs[a0].head;
            if !self.is_sink[j] || self.parent[j] == INVALID_IND {
                a0 = self.arcs[a0].next;
                continue;
            }
            // checking the origin of j
            let mut d = 0;
            loop {
                if self.ts[j] == self.time {
                    d += self.dist[j];
                    break;
                }
                let a = self.parent[j];
                d += 1;
                if a == TERMINAL {
                    self.ts[j] = self.time;
                    self.dist[j] = 1;
                    break;
                }
                if a == ORPHAN {
                    d = INVALID_IND;
                    break;
                }
                j = self.arcs[a].head;
            }
            if d != INVALID_IND {
                // j originates from the sink - done
                if d < d_min {
                    a0_min = a0;
                    d_min = d;
                }
                // set marks along the path
                j = self.arcs[a0].head;

                loop {
                    if self.ts[j] == self.time {
                        break;
                    }
                    self.ts[j] = self.time;
                    self.dist[j] = d;
                    d -= 1;
                    j = self.arcs[self.parent[j]].head;
                }
            }
            a0 = self.arcs[a0].next;
        }

        self.parent[i] = a0_min;
        if self.parent[i] != INVALID_IND {
            self.ts[i] = self.time;
            self.dist[i] = d_min + 1;
        } else {
            // no parent is found, process neighbors
            a0 = self.first_arc[i];
            while a0 != INVALID_IND {
                let j = self.arcs[a0].head;
                let a = self.parent[j];
                if self.is_sink[j] && a != INVALID_IND {
                    if self.arcs[a0].r_cap != 0.0 {
                        self.set_active(j);
                    }
                    if a != TERMINAL && a != ORPHAN && self.arcs[a].head == i {
                        self.set_orphan_rear(j); // add j to the end of the adoption list
                    }
                }
                a0 = self.arcs[a0].next;
            }
        }
    }

    fn augment(&mut self, middle_arc_id: usize) {
        let mut bottle_neck = self.arcs[middle_arc_id].r_cap;
        let middle_arc_sister_id = self.arcs[middle_arc_id].sister;
        let mut i = self.arcs[middle_arc_sister_id].head;
        let mut aid;
        loop {
            aid = self.parent[i];
            if aid == TERMINAL {
                break;
            }
            let sister = &self.arcs[self.arcs[aid].sister];
            if bottle_neck > sister.r_cap {
                bottle_neck = sister.r_cap;
            }
            i = self.arcs[aid].head;
        }

        if bottle_neck > self.tr_cap[i] {
            bottle_neck = self.tr_cap[i];
        }

        i = self.arcs[middle_arc_id].head;
        loop {
            aid = self.parent[i];
            if aid == TERMINAL {
                break;
            }
            let arc = &self.arcs[aid];
            if bottle_neck > arc.r_cap {
                bottle_neck = arc.r_cap;
            }
            i = arc.head;
        }

        if bottle_neck > -self.tr_cap[i] {
            bottle_neck = -self.tr_cap[i];
        }

        self.arcs[middle_arc_sister_id].r_cap += bottle_neck;
        self.arcs[middle_arc_id].r_cap -= bottle_neck;

        // the source tree
        i = self.arcs[middle_arc_sister_id].head;
        loop {
            aid = self.parent[i];
            if aid == TERMINAL {
                break;
            }

            self.arcs[aid].r_cap += bottle_neck;
            let sister_id = self.arcs[aid].sister;
            let sister = &mut self.arcs[sister_id];
            sister.r_cap -= bottle_neck;
            if sister.r_cap == 0.0 {
                self.set_orphan_front(i);
            }
            i = self.arcs[aid].head;
        }

        self.tr_cap[i] -= bottle_neck;
        if self.tr_cap[i] == 0.0 {
            self.set_orphan_front(i);
        }

        // the sink tree
        i = self.arcs[middle_arc_id].head;
        loop {
            aid = self.parent[i];
            if aid == TERMINAL {
                break;
            }
            let sister_id = self.arcs[aid].sister;
            self.arcs[sister_id].r_cap += bottle_neck;
            self.arcs[aid].r_cap -= bottle_neck;
            if self.arcs[aid].r_cap == 0.0 {
                self.set_orphan_front(i);
            }
            i = self.arcs[aid].head;
        }
        self.tr_cap[i] += bottle_neck;
        if self.tr_cap[i] == 0.0 {
            self.set_orphan_front(i);
        }
        self.flow += bottle_neck;
    }

    pub fn max_flow(&mut self) -> f64 {
        let mut current_node = INVALID_IND;
        loop {
            if self.time == 3218 {
                println!("here");
            }
            let mut i = current_node;
            if i != INVALID_IND {
                self.next[i] = INVALID_IND;
                if self.parent[i] == INVALID_IND {
                    i = INVALID_IND;
                }
            }
            if i == INVALID_IND {
                i = self.next_active();
                if i == INVALID_IND {
                    break;
                }
            }
            let mut aid;
            if !self.is_sink[i] {
                aid = self.first_arc[i];
                while aid != INVALID_IND {
                    if self.arcs[aid].r_cap != 0.0 {
                        let j = self.arcs[aid].head;
                        if self.parent[j] == INVALID_IND {
                            self.is_sink[j] = false;
                            self.parent[j] = self.arcs[aid].sister;
                            self.ts[j] = self.ts[i];
                            self.dist[j] = self.dist[i] + 1;
                            self.set_active(j);
                        } else if self.is_sink[j] {
                            break;
                        } else if self.ts[j] <= self.ts[i] && self.dist[j] > self.dist[i] {
                            self.parent[j] = self.arcs[aid].sister;
                            self.ts[j] = self.ts[i];
                            self.dist[j] = self.dist[i] + 1;
                        }
                    }
                    aid = self.arcs[aid].next;
                }
            } else {
                aid = self.first_arc[i];
                while aid != INVALID_IND {
                    let sister_id = self.arcs[aid].sister;
                    if self.arcs[sister_id].r_cap != 0.0 {
                        let j = self.arcs[aid].head;
                        if self.parent[j] == INVALID_IND {
                            self.is_sink[j] = true;
                            self.parent[j] = sister_id;
                            self.ts[j] = self.ts[i];
                            self.dist[j] = self.dist[i] + 1;
                            self.set_active(j);
                        } else if !self.is_sink[j] {
                            aid = sister_id;
                            break;
                        } else if self.ts[j] <= self.ts[i] && self.dist[j] > self.dist[i] {
                            self.parent[j] = sister_id;
                            self.ts[j] = self.ts[i];
                            self.dist[j] = self.dist[i] + 1;
                        }
                    }
                    aid = self.arcs[aid].next;
                }
            }
            self.time += 1;
            if aid != INVALID_IND {
                self.next[i] = i;
                current_node = i;
                self.augment(aid);
                i = self.orphan_first;
                while i != INVALID_IND {
                    let next_node = self.next_orphan[i];
                    self.next_orphan[i] = INVALID_IND;
                    i = self.orphan_first;
                    while i != INVALID_IND {
                        self.orphan_first = self.next_orphan[i];
                        if self.orphan_first == INVALID_IND {
                            self.orphan_last = INVALID_IND;
                        }
                        if self.is_sink[i] {
                            self.process_sink_orphan(i);
                        } else {
                            self.process_source_orphan(i);
                        }
                        i = self.orphan_first;
                    }
                    self.orphan_first = next_node;
                    i = self.orphan_first;
                }
            } else {
                current_node = INVALID_IND;
            }
        }
        let _count = self.count_sink();
        return self.flow;
    }

    fn count_sink(&self) -> usize {
        self.is_sink.iter().enumerate().filter_map(|(i, &v)| if v {Some(i)} else {None}).count()
    }
}
