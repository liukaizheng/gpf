use std::{
    alloc::Allocator,
    collections::{BinaryHeap, HashMap},
};

use bumpalo::{collections::CollectIn, Bump};
use itertools::Itertools;

use crate::{
    geometry::{Surf, Surface},
    mesh::EdgeId,
    point,
};

use super::TetSet;

struct EdgeAndLen {
    eid: EdgeId,
    len: f64,
}

impl PartialEq for EdgeAndLen {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.eid == other.eid && self.len == other.len
    }
}

impl Eq for EdgeAndLen {}

impl PartialOrd for EdgeAndLen {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.len.partial_cmp(&other.len)
    }
}

impl Ord for EdgeAndLen {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

struct SubdivisionData<'a> {
    surfaces: Vec<&'a Surf>,
    vals_and_grads: Vec<HashMap<usize, [f64; 4]>>,
    active_surfs: Vec<Vec<usize>>,
    queue: BinaryHeap<EdgeAndLen>,
}

pub(super) fn adaptive_subdivide(tets: &mut TetSet, surfaces: Vec<&Surf>) {
    let bump = Bump::new();
    let mut vals_and_grads = vec![HashMap::<usize, [f64; 4]>::new(); surfaces.len()];
    // TODO: parallelize by rayon
    for (i, &surf) in surfaces.iter().enumerate() {
        for (vid, p) in tets.points.chunks(3).enumerate() {
            vals_and_grads[i].insert(vid, surf.eval(p));
        }
    }
    let mut active_surfs = (0..tets.tets.len())
        .map(|_| Vec::from_iter(0..surfaces.len()))
        .collect_vec();
    let mut data = SubdivisionData {
        surfaces,
        vals_and_grads,
        active_surfs,
        queue: BinaryHeap::new(),
    };
    for tid in 0..tets.tets.len() {
        push_longest_edge(tid, tets, &mut data, &bump);
    }
}

fn push_longest_edge<A: Allocator + Copy>(
    tid: usize,
    tets: &mut TetSet,
    data: &mut SubdivisionData,
    alloc: A,
) {
}

fn is_subdivided<A: Allocator + Copy>(
    tid: usize,
    tets: &mut TetSet,
    data: &mut SubdivisionData,
    alloc: A,
) {
    let verts = tets.vertices_in(tid, alloc);
    let points = {
        let mut points = Vec::with_capacity_in(4, alloc);
        points.extend(verts.iter().map(|vid| point(&tets.points, vid.0)));
        points
    };

    for &sid /*surface id*/ in &data.active_surfs[tid] {
    }
}
