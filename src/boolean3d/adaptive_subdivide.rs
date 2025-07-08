use core::panic;
use std::{alloc::Allocator, collections::BinaryHeap};

use hashbrown::{HashMap, HashSet};

use bumpalo::Bump;
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    INVALID_IND,
    geometry::{BBox, Surf, Surface},
    math::{cross, cross_in, dot, interpolate, square_norm, sub_short},
    mesh::{EdgeId, Mesh},
    point, point_2, point_3,
    triangle::{convex_2, convex_3},
};

#[derive(Default)]
pub(crate) struct SurfaceEvaluation {
    pub(crate) sid: usize,
    pub(crate) evaluation: [[f64; 4]; 4],
}

#[derive(Clone, PartialEq, Eq)]
pub(crate) struct TetHandle {
    pub(crate) tid: usize,
    pub(crate) ver: usize,
}

impl Default for TetHandle {
    fn default() -> Self {
        Self {
            tid: INVALID_IND,
            ver: 12,
        }
    }
}

impl TetHandle {
    #[inline]
    pub fn new(tid: usize, ver: usize) -> Self {
        Self { tid, ver }
    }

    #[inline]
    pub(crate) fn sym_in_place(&mut self) -> &mut Self {
        self.ver = sym(self.ver);
        self
    }
}

#[inline]
fn next(ver: usize) -> usize {
    (ver + 4) % 12
}

#[inline]
fn prev(ver: usize) -> usize {
    (ver + 8) % 12
}

#[inline]
fn sym(ver: usize) -> usize {
    const SYM_TBL: [usize; 12] = [9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2];
    SYM_TBL[ver]
}

#[inline]
fn org(ver: usize) -> usize {
    const ORG_PIVOT: [usize; 12] = [3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0];
    ORG_PIVOT[ver]
}
#[inline]
fn dest(ver: usize) -> usize {
    const DEST_PIVOT: [usize; 12] = [2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1];
    DEST_PIVOT[ver]
}
#[inline]
fn apex(ver: usize) -> usize {
    const APEX_PIVOT: [usize; 12] = [1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2];
    APEX_PIVOT[ver]
}
#[inline]
fn oppo(ver: usize) -> usize {
    const OPPO_PIVOT: [usize; 12] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
    OPPO_PIVOT[ver]
}

#[inline(always)]
fn tet_edge_index(pa: usize, pb: usize) -> usize {
    pa + pb - (pa.min(pb) == 0) as usize
}
pub(crate) struct Tet {
    pub(crate) vertices: [usize; 4],
    pub(crate) neighbors: [TetHandle; 4],
    pub(crate) surface_evaluations: TinyVec<[SurfaceEvaluation; 3]>,
    pub(crate) edge_square_lengths: [f64; 6],
}

const VER_TO_EDGE: [usize; 12] = [5, 2, 0, 3, 3, 1, 2, 1, 4, 5, 4, 0];
impl Default for Tet {
    fn default() -> Self {
        Self {
            vertices: [INVALID_IND; 4],
            neighbors: [
                TetHandle::default(),
                TetHandle::default(),
                TetHandle::default(),
                TetHandle::default(),
            ],
            surface_evaluations: TinyVec::new(),
            edge_square_lengths: [0.0; 6],
        }
    }
}

impl Tet {
    #[inline]
    fn org(&self, hv: usize) -> usize {
        self.vertices[org(hv)]
    }
    #[inline]
    fn dest(&self, hv: usize) -> usize {
        self.vertices[dest(hv)]
    }
    #[inline]
    fn apex(&self, hv: usize) -> usize {
        self.vertices[apex(hv)]
    }
    #[inline]
    fn oppo(&self, hv: usize) -> usize {
        self.vertices[oppo(hv)]
    }

    #[inline]
    fn largest_edge(&self) -> (usize, f64) {
        let idx = self
            .edge_square_lengths
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        ([11, 5, 6, 3, 8, 9][idx], self.edge_square_lengths[idx])
    }

    #[inline]
    fn edge_square_len(&self, ver: usize) -> f64 {
        self.edge_square_lengths[VER_TO_EDGE[ver]]
    }

    fn subdividable(&mut self, points: &[f64], srf_datum: &[SurfaceData], sq_eps: f64) -> bool {
        if self.vertices[3] == INVALID_IND {
            return false;
        }

        let tet_points = self.vertices.map(|vid| point::<3>(points, vid));
        let tet_box = BBox::from_iter(tet_points);

        let mut contain_some_srf = false;
        self.surface_evaluations.retain(|eval| {
            let srf_data = &srf_datum[eval.sid];
            if tet_box.contains(&srf_data.bbox) {
                contain_some_srf = true;
                return true;
            }
            if !tet_box.intersects(&srf_data.bbox) {
                return false;
            }
            if srf_data.sub_bboxes.len() > 1 {
                if srf_data
                    .sub_bboxes
                    .iter()
                    .any(|bbox| tet_box.contains(bbox))
                {
                    contain_some_srf = true;
                    return true;
                }
                for bbox in &srf_data.sub_bboxes {
                    if tet_box.intersects(bbox) {
                        return true;
                    }
                }
                false
            } else {
                true
            }
        });
        if contain_some_srf {
            return true;
        }

        if self.surface_evaluations.len() < 1 {
            return false;
        }

        let trans_vmat = [
            sub_short::<3, _>(&tet_points[1], &tet_points[0]),
            sub_short::<3, _>(&tet_points[2], &tet_points[0]),
            sub_short::<3, _>(&tet_points[3], &tet_points[0]),
            sub_short::<3, _>(&tet_points[2], &tet_points[1]),
            sub_short::<3, _>(&tet_points[3], &tet_points[1]),
            sub_short::<3, _>(&tet_points[3], &tet_points[2]),
        ];

        let sq_det_vmat = {
            let d = det(&trans_vmat);
            d * d
        };
        let adj_vmat = [
            cross(&trans_vmat[1], &trans_vmat[2]),
            cross(&trans_vmat[2], &trans_vmat[0]),
            cross(&trans_vmat[0], &trans_vmat[1]),
        ];

        let n_surfaces = self.surface_evaluations.len();

        let mut interpolant_vec = Vec::with_capacity(n_surfaces);
        let mut interpolant_diff_vec = Vec::with_capacity(n_surfaces);
        let mut val_diff_vec = Vec::with_capacity(n_surfaces);
        for eval in self.surface_evaluations.iter() {
            let tet_vals_grads = &eval.evaluation;

            let mut vals = Vec::with_capacity(20);
            vals.extend(tet_vals_grads.iter().map(|vals_grads| vals_grads[0]));

            let v0 = tet_vals_grads[0][0];
            let g = &tet_vals_grads[0][1..];
            const S: f64 = 1.0 / 3.0;
            vals.extend([
                v0 + S * dot(g, &trans_vmat[0]),
                v0 + S * dot(g, &trans_vmat[1]),
                v0 + S * dot(g, &trans_vmat[2]),
            ]);

            let v1 = tet_vals_grads[1][0];
            let g = &tet_vals_grads[1][1..];
            vals.extend([
                v1 + S * dot(g, &trans_vmat[3]),
                v1 + S * dot(g, &trans_vmat[4]),
                v1 - S * dot(g, &trans_vmat[0]),
            ]);

            let v2 = tet_vals_grads[2][0];
            let g = &tet_vals_grads[2][1..];
            vals.extend([
                v2 + S * dot(g, &trans_vmat[5]),
                v2 - S * dot(g, &trans_vmat[1]),
                v2 - S * dot(g, &trans_vmat[3]),
            ]);

            let v3 = tet_vals_grads[3][0];
            let g = &tet_vals_grads[3][1..];
            vals.extend([
                v3 - S * dot(g, &trans_vmat[2]),
                v3 - S * dot(g, &trans_vmat[4]),
                v3 - S * dot(g, &trans_vmat[5]),
            ]);

            vals.push(
                ((vals[7] + vals[8] + vals[10] + vals[12] + vals[14] + vals[15]) * 1.5
                    - vals[1]
                    - vals[2]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[5] + vals[6] + vals[10] + vals[11] + vals[13] + vals[15]) * 1.5
                    - vals[0]
                    - vals[2]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[4] + vals[6] + vals[8] + vals[9] + vals[13] + vals[14]) * 1.5
                    - vals[0]
                    - vals[1]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[4] + vals[5] + vals[7] + vals[9] + vals[11] + vals[12]) * 1.5
                    - vals[0]
                    - vals[1]
                    - vals[2])
                    / 6.0,
            );

            let mut diffs = Vec::with_capacity(16);
            for i in 0..16 {
                let c = &C[i];
                diffs.push(v0 * c[0] + v1 * c[1] + v2 * c[2] + v3 * c[3] - vals[i + 4]);
            }

            let val_diff = [v1 - v0, v2 - v0, v3 - v0];
            let is_active = *vals
                .iter()
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap()
                > 0.0
                && *vals
                    .iter()
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap()
                    < 0.0;
            if is_active {
                if test_distance_1(&adj_vmat, val_diff, &diffs, sq_det_vmat, sq_eps) {
                    return true;
                }
                interpolant_vec.push(vals);
                interpolant_diff_vec.push(diffs);
                val_diff_vec.push(val_diff);
            }
        }

        if interpolant_vec.len() < 2 {
            return false;
        }

        let mut pair_set = HashSet::new();
        for ((i, v1), (j, v2)) in interpolant_vec.iter().enumerate().tuple_windows() {
            let mut points = Vec::with_capacity(42);
            points.extend(v1.iter().interleave(v2).map(|v| *v));

            if !contain_zero_2(points, std::alloc::Global) {
                continue;
            }

            pair_set.insert([i, j]);

            let h = [val_diff_vec[i], val_diff_vec[j]];
            let b = [
                interpolant_diff_vec[i].as_slice(),
                interpolant_diff_vec[j].as_slice(),
            ];
            if test_distance_2(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
                return true;
            }
        }

        for ((i, v1), (j, v2), (k, v3)) in interpolant_vec.iter().enumerate().tuple_windows() {
            if !pair_set.contains(&[i, j])
                || !pair_set.contains(&[i, k])
                || !pair_set.contains(&[j, k])
            {
                continue;
            }
            let mut points = Vec::with_capacity(63);
            points.extend(
                v1.iter()
                    .zip(v2)
                    .zip(v3)
                    .map(|((a, b), c)| [*a, *b, *c])
                    .flatten(),
            );

            if !contain_zero_3(points, std::alloc::Global) {
                continue;
            }

            let h = [val_diff_vec[i], val_diff_vec[j], val_diff_vec[k]];
            let b = [
                interpolant_diff_vec[i].as_slice(),
                interpolant_diff_vec[j].as_slice(),
                interpolant_diff_vec[k].as_slice(),
            ];
            if test_distance_3(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
                return true;
            }
        }
        return false;
    }

    fn subdivide(
        &mut self,
        edge: &TetHandle,
        new_tid: usize,
        points: &[f64],
        srf_datum: &[SurfaceData],
        length_map: &mut HashMap<usize, f64>,
        evaluation_map: &mut HashMap<usize, [f64; 4]>,
    ) -> Tet {
        let &TetHandle {tid, ver} = edge;
        let indices = [org(ver), dest(ver), apex(ver), oppo(ver)];
        let edge_map = [
            tet_edge_index(indices[0], indices[1]),
            tet_edge_index(indices[0], indices[2]),
            tet_edge_index(indices[0], indices[3]),
            tet_edge_index(indices[1], indices[2]),
            tet_edge_index(indices[1], indices[3]),
            tet_edge_index(indices[2], indices[3]),
        ];
        let vs: [usize; 5] = [
            self.vertices[indices[0]],
            self.vertices[indices[1]],
            self.vertices[indices[2]],
            self.vertices[indices[3]],
            points.len() / 3 - 1,
        ];

        let half_len = *length_map
            .entry(vs[0])
            .or_insert(self.edge_square_lengths[VER_TO_EDGE[ver]] * 0.25);
        let v2_v4_len = *length_map
            .entry(vs[2])
            .or_insert(square_norm(&sub_short::<3, _>(
                point::<3>(points, vs[2]),
                &points[(points.len() - 3)..],
            )));
        let v3_v4_len = *length_map
            .entry(vs[3])
            .or_insert(square_norm(&sub_short::<3, _>(
                point::<3>(points, vs[3]),
                &points[(points.len() - 3)..],
            )));

        let mut get_srf_eval = |indices: [usize; 4]| {
            TinyVec::from_iter(self.surface_evaluations.iter().map(|tet_eval| {
                let sid = tet_eval.sid;
                let evaluation = indices.map(|idx| {
                    if idx < 4 {
                        tet_eval.evaluation[idx]
                    } else {
                        let sid = tet_eval.sid;
                        *evaluation_map
                            .entry(sid)
                            .or_insert(srf_datum[sid].surf.eval(&points[(points.len() - 3)..]))
                    }
                });
                SurfaceEvaluation { sid, evaluation }
            }))
        };

        if vs[2] == INVALID_IND {
            let new_tet = {
                const VERT_INDICES: [usize; 4] = [1, 4, 3, 2];
                let vertices = VERT_INDICES.map(|i| vs[i]);
                let edge_square_lengths = [
                    half_len,
                    self.edge_square_lengths[4],
                    self.edge_square_lengths[3],
                    v3_v4_len,
                    v2_v4_len,
                    self.edge_square_lengths[5],
                ];
                Tet {
                    vertices,
                    edge_square_lengths,
                    surface_evaluations: get_srf_eval(VERT_INDICES),
                    neighbors: [
                        TetHandle::new(tid, 9),
                        self.neighbors[indices[0]].clone(),
                        TetHandle::default(),
                        TetHandle::default(),
                    ],
                }
            };
            {
                const VERT_INDICES: [usize; 4] = [4, 0, 3, 2];
                self.vertices = VERT_INDICES.map(|i| vs[i]);
                self.edge_square_lengths = [
                    half_len,
                    v3_v4_len,
                    v2_v4_len,
                    self.edge_square_lengths[2],
                    self.edge_square_lengths[1],
                    self.edge_square_lengths[5],
                ];
                self.surface_evaluations = get_srf_eval(VERT_INDICES);
                self.neighbors = [
                    self.neighbors[indices[1]].clone(),
                    TetHandle::new(new_tid, 8),
                    TetHandle::default(),
                    TetHandle::default(),
                ];
            }
            new_tet
        } else {
            let new_tet = {
                const VERT_INDICES: [usize; 4] = [4, 1, 2, 3];
                let vertices = VERT_INDICES.map(|i| vs[i]);
                let edge_square_lengths = [
                    half_len,
                    v2_v4_len,
                    v3_v4_len,
                    self.edge_square_lengths[3],
                    self.edge_square_lengths[4],
                    self.edge_square_lengths[5],
                ];
                Tet {
                    vertices,
                    edge_square_lengths,
                    surface_evaluations: get_srf_eval(VERT_INDICES),
                    neighbors: [
                        self.neighbors[indices[0]].clone(),
                        TetHandle::new(tid, 8),
                        Default::default(),
                        Default::default(),
                    ]
                }
            };
            {
                const VERT_INDICES: [usize; 4] = [0, 4, 2, 3];
                self.vertices = VERT_INDICES.map(|i| vs[i]);
                self.edge_square_lengths = [
                    half_len,
                    self.edge_square_lengths[1],
                    self.edge_square_lengths[2],
                    v2_v4_len,
                    v3_v4_len,
                    self.edge_square_lengths[5],
                ];
                self.neighbors = [
                    TetHandle::new(new_tid, 9),
                    self.neighbors[indices[1]].clone(),
                    TetHandle::default(),
                    TetHandle::default(),
                ];
            }
            new_tet
        }
    }
}

pub(crate) struct TetComplex {
    pub(crate) points: Vec<f64>,
    pub(crate) tets: Vec<Tet>,
}

impl TetComplex {
    pub(crate) fn spin_edges(&self, first_handle: TetHandle) -> impl Iterator<Item = TetHandle> {
        let mut next_handle = first_handle.clone();
        let mut is_first = true;
        std::iter::from_fn(move || {
            if !is_first && next_handle == first_handle {
                None
            } else {
                is_first = false;
                let curr_handle = next_handle.clone();
                self.sym_face_in_place(next_handle.sym_in_place());
                Some(curr_handle)
            }
        })
    }

    pub(crate) fn sym_face_in_place(&self, handle: &mut TetHandle) {
        const FSYM_TBL: [[usize; 12]; 12] = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7],
            [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7],
            [4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3],
            [4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3],
            [4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3],
            [4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3],
        ];
        let nei = &self.tets[handle.tid].neighbors[handle.tid & 3];
        handle.tid = nei.tid;
        handle.ver = FSYM_TBL[handle.ver][nei.ver];
    }
}

use super::TetSet;

#[derive(Debug, Clone)]
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

pub(crate) struct SurfaceData<'a> {
    pub(crate) surf: &'a Surf,
    pub(crate) bbox: BBox,
    pub(crate) sub_bboxes: TinyVec<[BBox; 1]>,
}

struct SubdivisionData<'a> {
    surface_datum: Vec<SurfaceData<'a>>,
    vals_and_grads: Vec<Vec<[f64; 4]>>,
    queue: BinaryHeap<EdgeAndLen>,
}

struct TetEdgeAndLen {
    handle: TetHandle,
    len: f64,
}

impl PartialEq for TetEdgeAndLen {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.handle == other.handle && self.len == other.len
    }
}

impl Eq for TetEdgeAndLen {}

impl PartialOrd for TetEdgeAndLen {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.len.partial_cmp(&other.len)
    }
}

impl Ord for TetEdgeAndLen {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

pub(super) fn adaptive_subdivide1<'a>(
    tets: &mut TetComplex,
    srf_datum: Vec<SurfaceData<'a>>,
    sq_eps: f64,
) {
    let mut pq = BinaryHeap::new();
    for (tid, tet) in &mut tets.tets[..6].iter_mut().enumerate() {
        if tet.subdividable(&tets.points, &srf_datum, sq_eps) {
            let (ver, len) = tet.largest_edge();
            pq.push(TetEdgeAndLen {
                handle: TetHandle { tid, ver },
                len,
            });
        }
    }

    let mut spin_edges: Vec<TetHandle> = Vec::new();
    let mut split_tets: Vec<[TetHandle; 2]> = Vec::new();
    let mut length_map = HashMap::new();
    let mut evaluation_map = HashMap::with_capacity(srf_datum.len());

    loop {
        if pq.is_empty() {
            break;
        }

        let TetEdgeAndLen { handle: curr, len } = pq.pop().unwrap();
        let tet = &mut tets.tets[curr.tid];
        if len != tet.edge_square_len(curr.ver) {
            continue;
        }

        let va = tet.org(curr.ver);
        let vb = tet.dest(curr.ver);
        tets.points.extend_from_slice(&interpolate::<3>(
            point::<3>(&tets.points, va),
            point::<3>(&tets.points, vb),
            0.5,
        ));
        spin_edges.extend(tets.spin_edges(curr));
        split_tets.extend(spin_edges.iter().map(|handle| {
            let new_tid = tets.tets.len();
            let new_tet = tets.tets[handle.tid].subdivide(
                handle,
                new_tid,
                &tets.points,
                &srf_datum,
                &mut length_map,
                &mut evaluation_map,
            );
            tets.tets.push(new_tet);
            [TetHandle::default(), TetHandle::default()]
        }));
        spin_edges.clear();
        length_map.clear();
        evaluation_map.clear();
    }
}

pub(super) fn adaptive_subdivide<'a>(
    tets: &mut TetSet,
    surface_datum: Vec<SurfaceData<'a>>,
    sq_eps: f64,
) -> Vec<Vec<f64>> {
    let mut vals_and_grads = vec![Vec::with_capacity(tets.mesh.n_vertices()); surface_datum.len()];
    for p in tets.points.chunks(3) {
        for (i, surf) in surface_datum.iter().map(|d| &d.surf).enumerate() {
            vals_and_grads[i].push(surf.eval(p));
        }
    }

    let mut data = SubdivisionData {
        surface_datum,
        vals_and_grads,
        queue: BinaryHeap::new(),
    };

    adaptive_subdivide_impl(tets, &mut data, sq_eps);

    data.vals_and_grads
        .into_iter()
        .map(|v| v.into_iter().map(|vg| vg[0]).collect_vec())
        .collect_vec()
}

fn adaptive_subdivide_impl(tets: &mut TetSet, data: &mut SubdivisionData, sq_eps: f64) {
    let mut check_bump = Bump::new();
    for tid in 0..tets.tet_faces.len() {
        check_bump.reset();
        push_longest_edge(tid, tets, data, sq_eps, &check_bump);
    }

    let mut split_bump = Bump::new();
    while !data.queue.is_empty() {
        let EdgeAndLen { eid, len } = data.queue.pop().unwrap();
        if tets.square_edge_lengths[eid] != len {
            // edge changed
            continue;
        }
        split_bump.reset();
        let (new_vert, tet_pairs) = tets.split_edge(eid, &split_bump);
        let p = point_3(&tets.points, new_vert.0);

        for (sid, surf) in data.surface_datum.iter().map(|d| &d.surf).enumerate() {
            data.vals_and_grads[sid].push(surf.eval(p));
        }

        for tid in tet_pairs.into_iter().flatten() {
            check_bump.reset();
            push_longest_edge(tid, tets, data, sq_eps, &check_bump);
        }
    }
}

fn push_longest_edge(
    tid: usize,
    tets: &mut TetSet,
    data: &mut SubdivisionData,
    sq_eps: f64,
    bump: &Bump,
) {
    if tets.tet_edges[tid]
        .iter()
        .all(|&eid| tets.square_edge_lengths[eid] < sq_eps)
    {
        return;
    }
    if subdividable(tid, tets, data, sq_eps, bump) {
        let longest_eid = *tets.tet_edges[tid]
            .iter()
            .max_by(|&&ea, &&eb| {
                tets.square_edge_lengths[ea]
                    .partial_cmp(&tets.square_edge_lengths[eb])
                    .unwrap()
            })
            .unwrap();
        data.queue.push(EdgeAndLen {
            eid: longest_eid,
            len: tets.square_edge_lengths[longest_eid],
        });
    }
}

const C: [[f64; 4]; 16] = [
    [2.0 / 3.0, 1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 1.0 / 3.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0, 1.0 / 3.0],
    [0.0, 2.0 / 3.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0],
    [0.0, 0.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 2.0 / 3.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 0.0, 2.0 / 3.0],
    [0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
];

fn subdividable<A: Allocator + Copy>(
    tid: usize,
    tets: &mut TetSet,
    data: &SubdivisionData,
    sq_eps: f64,
    alloc: A,
) -> bool {
    let verts = tets.tet_vertices[tid];
    let tet_points = verts.map(|vid| point_3(&tets.points, vid.0));
    let tet_box = BBox::from_iter(tet_points);

    let mut contain_some_srf = false;

    let tet_surfaces = &mut tets.surf_indices[tid];
    tet_surfaces.retain(|&sid| {
        let srf_data = &data.surface_datum[sid];
        if tet_box.contains(&srf_data.bbox) {
            contain_some_srf = true;
            return true;
        }
        if !tet_box.intersects(&srf_data.bbox) {
            return false;
        }
        if srf_data.sub_bboxes.len() > 1 {
            if srf_data
                .sub_bboxes
                .iter()
                .any(|bbox| tet_box.contains(bbox))
            {
                contain_some_srf = true;
                return true;
            }
            for bbox in &srf_data.sub_bboxes {
                if tet_box.intersects(bbox) {
                    return true;
                }
            }
            false
        } else {
            true
        }
    });

    if contain_some_srf {
        return true;
    }

    if tet_surfaces.len() < 1 {
        return false;
    }

    let trans_vmat = [
        sub_short::<3, _>(tet_points[1], tet_points[0]),
        sub_short::<3, _>(tet_points[2], tet_points[0]),
        sub_short::<3, _>(tet_points[3], tet_points[0]),
        sub_short::<3, _>(tet_points[2], tet_points[1]),
        sub_short::<3, _>(tet_points[3], tet_points[1]),
        sub_short::<3, _>(tet_points[3], tet_points[2]),
    ];

    let sq_det_vmat = {
        let d = det(&trans_vmat);
        d * d
    };

    let adj_vmat = [
        cross(&trans_vmat[1], &trans_vmat[2]),
        cross(&trans_vmat[2], &trans_vmat[0]),
        cross(&trans_vmat[0], &trans_vmat[1]),
    ];

    let n_surfaces = tet_surfaces.len();

    let mut interpolant_vec = Vec::with_capacity_in(n_surfaces, alloc);
    let mut interpolant_diff_vec = Vec::with_capacity_in(n_surfaces, alloc);
    let mut val_diff_vec = Vec::with_capacity_in(n_surfaces, alloc);
    for &sid in tet_surfaces.iter() {
        let tet_vals_grads = verts.map(|vid| &data.vals_and_grads[sid][vid]);

        let mut vals = Vec::with_capacity_in(20, alloc);
        vals.extend(tet_vals_grads.iter().map(|vals_grads| vals_grads[0]));

        let v0 = tet_vals_grads[0][0];
        let g = &tet_vals_grads[0][1..];
        const S: f64 = 1.0 / 3.0;
        vals.extend([
            v0 + S * dot(g, &trans_vmat[0]),
            v0 + S * dot(g, &trans_vmat[1]),
            v0 + S * dot(g, &trans_vmat[2]),
        ]);

        let v1 = tet_vals_grads[1][0];
        let g = &tet_vals_grads[1][1..];
        vals.extend([
            v1 + S * dot(g, &trans_vmat[3]),
            v1 + S * dot(g, &trans_vmat[4]),
            v1 - S * dot(g, &trans_vmat[0]),
        ]);

        let v2 = tet_vals_grads[2][0];
        let g = &tet_vals_grads[2][1..];
        vals.extend([
            v2 + S * dot(g, &trans_vmat[5]),
            v2 - S * dot(g, &trans_vmat[1]),
            v2 - S * dot(g, &trans_vmat[3]),
        ]);

        let v3 = tet_vals_grads[3][0];
        let g = &tet_vals_grads[3][1..];
        vals.extend([
            v3 - S * dot(g, &trans_vmat[2]),
            v3 - S * dot(g, &trans_vmat[4]),
            v3 - S * dot(g, &trans_vmat[5]),
        ]);

        vals.push(
            ((vals[7] + vals[8] + vals[10] + vals[12] + vals[14] + vals[15]) * 1.5
                - vals[1]
                - vals[2]
                - vals[3])
                / 6.0,
        );
        vals.push(
            ((vals[5] + vals[6] + vals[10] + vals[11] + vals[13] + vals[15]) * 1.5
                - vals[0]
                - vals[2]
                - vals[3])
                / 6.0,
        );
        vals.push(
            ((vals[4] + vals[6] + vals[8] + vals[9] + vals[13] + vals[14]) * 1.5
                - vals[0]
                - vals[1]
                - vals[3])
                / 6.0,
        );
        vals.push(
            ((vals[4] + vals[5] + vals[7] + vals[9] + vals[11] + vals[12]) * 1.5
                - vals[0]
                - vals[1]
                - vals[2])
                / 6.0,
        );

        let mut diffs = Vec::with_capacity_in(16, alloc);
        for i in 0..16 {
            let c = &C[i];
            diffs.push(v0 * c[0] + v1 * c[1] + v2 * c[2] + v3 * c[3] - vals[i + 4]);
        }

        let val_diff = [v1 - v0, v2 - v0, v3 - v0];
        let is_active = *vals
            .iter()
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
            > 0.0
            && *vals
                .iter()
                .min_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap()
                < 0.0;
        if is_active {
            if test_distance_1(&adj_vmat, val_diff, &diffs, sq_det_vmat, sq_eps) {
                return true;
            }
            interpolant_vec.push(vals);
            interpolant_diff_vec.push(diffs);
            val_diff_vec.push(val_diff);
        }
    }

    if interpolant_vec.len() < 2 {
        return false;
    }

    let mut pair_set = HashSet::new_in(alloc);
    for ((i, v1), (j, v2)) in interpolant_vec.iter().enumerate().tuple_windows() {
        let mut points = Vec::with_capacity_in(42, alloc);
        points.extend(v1.iter().interleave(v2).map(|v| *v));

        if !contain_zero_2(points, alloc) {
            continue;
        }

        pair_set.insert([i, j]);

        let h = [val_diff_vec[i], val_diff_vec[j]];
        let b = [
            interpolant_diff_vec[i].as_slice(),
            interpolant_diff_vec[j].as_slice(),
        ];
        if test_distance_2(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
            return true;
        }
    }

    for ((i, v1), (j, v2), (k, v3)) in interpolant_vec.iter().enumerate().tuple_windows() {
        if !pair_set.contains(&[i, j]) || !pair_set.contains(&[i, k]) || !pair_set.contains(&[j, k])
        {
            continue;
        }
        let mut points = Vec::with_capacity_in(63, alloc);
        points.extend(
            v1.iter()
                .zip(v2)
                .zip(v3)
                .map(|((a, b), c)| [*a, *b, *c])
                .flatten(),
        );

        if !contain_zero_3(points, alloc) {
            continue;
        }

        let h = [val_diff_vec[i], val_diff_vec[j], val_diff_vec[k]];
        let b = [
            interpolant_diff_vec[i].as_slice(),
            interpolant_diff_vec[j].as_slice(),
            interpolant_diff_vec[k].as_slice(),
        ];
        if test_distance_3(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
            return true;
        }
    }
    return false;
}

fn transpose_adjacent_mat<const N: usize>(mat: &[[f64; N]]) -> [[f64; N]; N] {
    let mut vec = [[0.0; N]; N];
    if N == 2 {
        vec[0][0] = mat[1][1];
        vec[0][1] = -mat[0][1];
        vec[1][0] = -mat[1][0];
        vec[1][1] = mat[0][0];
    } else if N == 3 {
        cross_in(&mat[1], &mat[2], &mut vec[0]);
        cross_in(&mat[2], &mat[0], &mut vec[1]);
        cross_in(&mat[0], &mat[1], &mut vec[2]);
    } else {
        panic!("not implemented");
    }
    vec
}

fn det<const N: usize>(mat: &[[f64; N]]) -> f64 {
    if N == 2 {
        mat[0][0] * mat[1][1] - mat[0][1] * mat[0][1]
    } else if N == 3 {
        mat[0][0] * mat[1][1] * mat[2][2]
            + mat[0][1] * mat[1][2] * mat[2][0]
            + mat[0][2] * mat[1][0] * mat[2][1]
            - mat[0][2] * mat[1][1] * mat[2][0]
            - mat[0][1] * mat[1][0] * mat[2][2]
            - mat[0][0] * mat[1][2] * mat[2][1]
    } else {
        panic!("not implemented");
    }
}

fn contain_zero_2<A: Allocator + Copy>(mut points: Vec<f64, A>, alloc: A) -> bool {
    points.extend([0.0, 0.0]);
    let zero_vid = points.len() >> 1;
    let hull = convex_2(&points, alloc);
    for vid in hull {
        let p = point_2(&points, vid);
        if vid == zero_vid || (p[0] == 0.0 && p[1] == 0.0) {
            return false;
        }
    }
    true
}

fn contain_zero_3<A: Allocator + Copy>(mut points: Vec<f64, A>, alloc: A) -> bool {
    let zero_vid = points.len() / 3;
    points.extend([0.0, 0.0, 0.0]);
    match convex_3(&points, true, alloc) {
        crate::triangle::Convex3Result::Dim3(hull) => {
            for vid in hull {
                let p = point_3(&points, vid);
                if vid == zero_vid || (p[0] == 0.0 && p[1] == 0.0 && p[2] == 0.0) {
                    return false;
                }
            }
            true
        }
        _ => false,
    }
}

fn test_distance_1(adj_v: &[[f64; 3]], h: [f64; 3], b: &[f64], sq_det_v: f64, sq_eps: f64) -> bool {
    // w: (M, 3)
    let mut w = [0.0f64; 3];
    for i in 0..3 {
        w[i] = h[0] * adj_v[0][i] + h[1] * adj_v[1][i] + h[2] * adj_v[2][i];
    }
    let w2 = square_norm(&w);
    let max_b = b
        .iter()
        .map(|s| s.abs())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let b2 = max_b * max_b;
    return b2 * sq_det_v > w2 * sq_eps;
}

fn test_distance_2(
    adj_v: &[[f64; 3]],
    h: &[[f64; 3]],
    b: [&[f64]; 2],
    sq_det_v: f64,
    sq_eps: f64,
) -> bool {
    // w: (M, 3)
    let mut w = [[0.0f64; 3]; 2];
    for i in 0..2 {
        for j in 0..3 {
            w[i][j] = h[i][0] * adj_v[0][j] + h[i][1] * adj_v[1][j] + h[i][2] * adj_v[2][j];
        }
    }
    // u = w * w^T with shape(M, M)
    let mut u = [[0.0; 2]; 2];

    for i in 0..2 {
        for j in 0..2 {
            u[i][j] = w[i][0] * w[j][0] + w[i][1] * w[j][1] + w[i][2] * w[j][2];
        }
    }

    let det_u = det(&u);

    // adj_u: (M, M)
    let trans_adj_u = transpose_adjacent_mat(&u);
    // wu = w^T x adj_u with shape (3, M)
    let mut wu = [[0.0; 2]; 3];
    for i in 0..3 {
        for j in 0..2 {
            for k in 0..2 {
                wu[i][j] += w[k][i] * trans_adj_u[j][k];
            }
        }
    }
    let r2 = (0..b[0].len())
        .map(|l| {
            let mut d = [0.0; 3];
            for i in 0..3 {
                for j in 0..2 {
                    d[i] += wu[i][j] * b[j][l];
                }
            }
            square_norm(&d)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let sq_det_u = det_u * det_u;

    return r2 * sq_det_v > sq_det_u * sq_eps;
}

fn test_distance_3(
    adj_v: &[[f64; 3]],
    h: &[[f64; 3]],
    b: [&[f64]; 3],
    sq_det_v: f64,
    sq_eps: f64,
) -> bool {
    // w: (M, 3)
    let mut w = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            w[i][j] = h[i][0] * adj_v[0][j] + h[i][1] * adj_v[1][j] + h[i][2] * adj_v[2][j];
        }
    }
    let det_w = det(&w);
    let trans_adj_w = transpose_adjacent_mat(&w);
    let r2 = (0..b[0].len())
        .map(|l| {
            let mut d = [0.0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    d[i] += trans_adj_w[i][j] * b[j][l];
                }
            }
            square_norm(&d)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    return r2 * sq_det_v > det_w * det_w * sq_eps;
}
