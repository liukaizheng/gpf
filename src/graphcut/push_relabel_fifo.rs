use std::{
    collections::VecDeque,
    ops::{AddAssign, SubAssign},
};

use hashbrown::HashMap;

use super::MaxFlow;
#[derive(Clone)]
struct Vertex<T> {
    label: usize,
    excess: T,
}

#[derive(Clone)]
struct Arc<T> {
    dist: usize,
    /// reverse arc index
    reverse: usize,
    capacity: T,
}

impl<T: Copy> Arc<T> {
    fn new(dist: usize, reverse: usize, capacity: T) -> Self {
        Self {
            dist,
            reverse,
            capacity,
        }
    }
    fn reverse_capacity(&self, container: &[Vec<Self>]) -> T {
        container[self.dist][self.reverse].capacity
    }
}

const ALPHA: usize = 6;
const BETA: usize = 12;

pub struct PushRelabelFifo<T> {
    residual_network: Vec<Vec<Arc<T>>>,
    vertices: Vec<Vertex<T>>,
    queue: VecDeque<usize>,
    dist_queue: VecDeque<(usize, usize)>,
    source: usize,
    sink: usize,
    n_pushes: usize,
    n_relabels: usize,
    relabel_progress: usize,
    not_reached: usize,
    n_global_relabels: usize,
    relabel_threshold: usize,
}

impl<T> PushRelabelFifo<T>
where
    T: Copy + PartialOrd + From<i32> + AddAssign + SubAssign,
{
    fn init(&mut self) {
        for i in 0..self.residual_network[self.source].len() {
            let arc = self.residual_network[self.source][i].clone();
            if arc.capacity > T::from(0) {
                self.vertices[arc.dist].excess = arc.capacity;
                self.residual_network[arc.dist][arc.reverse].capacity += arc.capacity;
                self.residual_network[self.source][i].capacity = T::from(0);
                // self.queue.push_back(arc.dist);
                self.n_pushes += 1;
            }
        }
        let sum = self.residual_network.iter().map(|v| v.len()).sum::<usize>();
        self.relabel_threshold = sum / 2 + self.residual_network.len() * ALPHA;
    }

    fn global_relabel(&mut self) {
        self.n_global_relabels += 1;
        for v in self.vertices.iter_mut() {
            v.label = self.not_reached;
        }

        self.queue.clear();
        self.dist_queue.clear();
        self.dist_queue.push_back((self.sink, 0));
        self.vertices[self.sink].label = 0;

        while !self.dist_queue.is_empty() {
            let (v, d) = self.dist_queue.pop_front().unwrap();
            for arc in &self.residual_network[v] {
                if arc.reverse_capacity(&self.residual_network) > T::from(0)
                    && self.vertices[arc.dist].label == self.not_reached
                {
                    let vert = &mut self.vertices[arc.dist];
                    vert.label = d + 1;
                    self.dist_queue.push_back((arc.dist, d + 1));
                    if vert.excess > T::from(0) {
                        self.queue.push_back(arc.dist);
                    }
                }
            }
        }
    }

    fn discharge(&mut self, v: usize) {
        let mut label = self.vertices[v].label;
        while label < self.not_reached {
            if self.push(v, label) {
                return;
            }
            label = self.relabel(v);
        }
    }

    fn push(&mut self, vid: usize, label: usize) -> bool {
        let target_label = label - 1;
        for i in 0..self.residual_network[vid].len() {
            let arc = self.residual_network[vid][i].clone();
            if arc.capacity > T::from(0) && self.vertices[arc.dist].label == target_label {
                self.n_pushes += 1;
                let delta = if arc.capacity.partial_cmp(&self.vertices[vid].excess)
                    == Some(std::cmp::Ordering::Less)
                {
                    arc.capacity
                } else {
                    self.vertices[vid].excess
                };
                if self.vertices[arc.dist].excess == T::from(0) && arc.dist != self.sink {
                    self.queue.push_back(arc.dist);
                }

                self.vertices[vid].excess -= delta;
                self.vertices[arc.dist].excess += delta;
                self.residual_network[vid][i].capacity -= delta;
                self.residual_network[arc.dist][arc.reverse].capacity += delta;
                if self.vertices[vid].excess == T::from(0) {
                    return true;
                }
            }
        }
        false
    }

    fn relabel(&mut self, vid: usize) -> usize {
        self.n_relabels += 1;
        self.relabel_progress += BETA;
        let new_label = self.calc_new_label(vid);
        self.vertices[vid].label = new_label;
        new_label
    }

    fn calc_new_label(&self, vid: usize) -> usize {
        let mut min_label = self.not_reached - 1;
        for arc in &self.residual_network[vid] {
            if arc.capacity > T::from(0) {
                min_label = min_label.min(self.vertices[arc.dist].label);
            }
        }
        min_label + 1
    }
}

impl<U, T> From<(U, usize)> for PushRelabelFifo<T>
where
    U: AsRef<[(usize, usize, T)]>,
    T: Copy + PartialOrd + From<i32> + AddAssign + SubAssign,
{
    fn from(data: (U, usize)) -> Self {
        let (arcs, n_nodes) = data;
        let mut edges: Vec<(usize, T)> = Vec::with_capacity(arcs.as_ref().len() * 2);
        let mut edge_map: Vec<HashMap<usize, usize>> = vec![HashMap::new(); n_nodes];
        for &(from, to, capacity) in arcs.as_ref() {
            if let Some(&edge_id) = edge_map[from].get(&to) {
                edges[edge_id].1 += capacity;
            } else {
                if let Some(&sister_edge_id) = edge_map[to].get(&from) {
                    edges[sister_edge_id + 1].1 += capacity;
                } else {
                    edge_map[from].insert(to, edges.len());
                    edges.push((to, capacity));
                    edges.push((from, T::from(0)));
                }
            }
        }

        let mut residual_network = vec![Vec::new(); n_nodes];
        for chunk in edges.chunks(2) {
            let ea = chunk[0];
            let eb = chunk[1];
            let ea_idx = residual_network[eb.0].len();
            let eb_idx = residual_network[ea.0].len();
            residual_network[eb.0].push(Arc::new(ea.0, eb_idx, ea.1));
            residual_network[ea.0].push(Arc::new(eb.0, ea_idx, eb.1));
        }
        let mut graph = Self {
            residual_network,
            vertices: vec![
                Vertex {
                    label: 0,
                    excess: T::from(0)
                };
                n_nodes
            ],
            queue: VecDeque::new(),
            dist_queue: VecDeque::new(),
            source: 0,
            sink: n_nodes - 1,
            n_pushes: 0,
            n_relabels: 0,
            relabel_progress: 0,
            not_reached: n_nodes,
            n_global_relabels: 0,
            relabel_threshold: 0,
        };
        graph.init();
        graph
    }
}

impl<T> MaxFlow for PushRelabelFifo<T>
where
    T: Copy + PartialOrd + From<i32> + AddAssign + SubAssign,
{
    fn find_max_flow(&mut self) {
        self.global_relabel();
        while !self.queue.is_empty() {
            let v = self.queue.pop_front().unwrap();
            self.discharge(v);
            if self.relabel_progress >= self.relabel_threshold * 2 {
                self.relabel_progress = 0;
                self.global_relabel();
            }
        }
    }

    #[inline(always)]
    fn is_sink(&self, i: usize) -> bool {
        self.vertices[i].label != self.not_reached
    }
}
