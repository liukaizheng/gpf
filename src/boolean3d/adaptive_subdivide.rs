use std::{
    alloc::Allocator,
    collections::{BinaryHeap, HashMap},
};

use bumpalo::Bump;
use itertools::Itertools;

use crate::{
    geometry::{Surf, Surface},
    math::{cross, dot, sub_short},
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

fn subdividable(tid: usize, tets: &mut TetSet, data: &mut SubdivisionData, bump: &Bump) {
    let verts = tets.vertices_in(tid, bump);
    let points = bumpalo::collections::Vec::from_iter_in(
        verts.iter().map(|vid| point(&tets.points, vid.0)),
        bump,
    );
    let vmat = bumpalo::vec![in &bump;
        sub_short::<3, _>(points[1], points[0]),
        sub_short::<3, _>(points[2], points[0]),
        sub_short::<3, _>(points[3], points[0]),
        sub_short::<3, _>(points[2], points[1]),
        sub_short::<3, _>(points[3], points[1]),
        sub_short::<3, _>(points[3], points[2]),
    ];

    let adj_vmat = bumpalo::vec![in &bump;
        cross(&vmat[1], &vmat[2]),
        cross(&vmat[2], &vmat[0]),
        cross(&vmat[0], &vmat[1]),
    ];

    let mut interpolant_vec = Vec::with_capacity_in(data.active_surfs.len(), bump);
    for &sid /*surface id*/ in &data.active_surfs[tid] {
        let tet_vals_grads = bumpalo::collections::Vec::from_iter_in(
            verts.iter().map(|&vid| data.vals_and_grads[sid].get(&vid.0).unwrap()), bump);
        let mut vals = Vec::with_capacity_in(20, bump);
        
        vals.extend(tet_vals_grads.iter().map(|vals_grads| vals_grads[0]));
        let v = tet_vals_grads[0][0];
        let g = &tet_vals_grads[0][1..];
            
        vals.extend([v + dot(g, &vmat[0]), v + dot(g, &vmat[1]), v + dot(g, &vmat[2])]);
        let v = tet_vals_grads[1][0];
        let g = &tet_vals_grads[1][1..];
        vals.extend([v + dot(g, &vmat[3]), v + dot(g, &vmat[4]), v - dot(g, &vmat[0])]);
        
        let v = tet_vals_grads[2][0];
        let g = &tet_vals_grads[2][1..];
        vals.extend([v + dot(g, &vmat[5]), v - dot(g, &vmat[1]), v - dot(g, &vmat[3])]);
        
        
        let v = tet_vals_grads[3][0];
        let g = &tet_vals_grads[3][1..];
        vals.extend([v - dot(g, &vmat[2]), v - dot(g, &vmat[4]), v - dot(g, &vmat[5])]);
        
        interpolant_vec.push(vals);
    }
}
fn test_distance<const M: usize>(v: &[[f64; 3]], h: &[[f64; M]]) -> bool {
    false
}
