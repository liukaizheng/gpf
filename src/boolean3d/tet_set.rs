use std::{alloc::Allocator, collections::BinaryHeap, ops::Deref};

use bumpalo::Bump;
use hashbrown::{HashMap, HashSet};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    INVALID_IND, decode_index,
    geometry::{BBox, Surf, Surface},
    math::{cross, cross_in, dot, square_norm, sub_short},
    mesh::{EdgeId, ElementId, FaceId, Halfedge, HalfedgeId, Mesh, SurfaceMesh, VertexId},
    point, point_3,
    triangle::{convex_2, convex_3},
    twin_index,
};

#[derive(Default)]
pub(crate) struct SurfaceEvaluation {
    pub(crate) sid: usize,
    pub(crate) evaluation: [[f64; 4]; 4],
}

pub(crate) struct SurfaceData<'a> {
    pub(crate) surf: &'a Surf,
    pub(crate) bbox: BBox,
    pub(crate) sub_bboxes: TinyVec<[BBox; 1]>,
}

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

pub(crate) struct Tet {
    pub(crate) vertices: [VertexId; 4],
    pub(crate) edges: [EdgeId; 6],
    pub(crate) faces: [FaceId; 4],
    pub(crate) surface_evaluations: TinyVec<[SurfaceEvaluation; 3]>,
}

impl Tet {
    #[inline]
    fn face_from_vertex(&self, vid: VertexId) -> FaceId {
        let idx = self.vertices.iter().position(|&v| v == vid).unwrap();
        self.faces[idx]
    }

    #[inline]
    fn face_from_edge(&self, eid: EdgeId, va: VertexId) -> [FaceId; 2] {
        const EDGE_VERTICES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
        let idx = self.edges.iter().position(|&e| e == eid).unwrap();
        let [i1, i2] = EDGE_VERTICES[idx];
        if self.vertices[i1] == va {
            [self.faces[i2], self.faces[i1]]
        } else {
            [self.faces[i1], self.faces[i2]]
        }
    }

    fn subdividable(
        &mut self,
        points: &[f64],
        edge_square_lengths: &[f64],
        srf_datum: &[SurfaceData],
        sq_eps: f64,
    ) -> bool {
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
        if !self.vertices[3].valid()
            || self
                .edges
                .into_iter()
                .all(|eid| edge_square_lengths[eid] < sq_eps)
        {
            return false;
        }

        let tet_points = self.vertices.map(|vid| point::<3>(points, *vid));
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
}

pub(crate) struct TetSet {
    pub(crate) points: Vec<f64>,
    pub(crate) mesh: SurfaceMesh<std::alloc::Global>,
    pub(crate) tets: Vec<Tet>,
    pub(crate) face_tets: Vec<[usize; 2]>,
    pub(crate) square_edge_lengths: Vec<f64>,
}

#[inline]
pub(crate) fn tet_face_reversed(face_tets: &[usize; 2], tid: usize) -> bool {
    debug_assert!(face_tets[0] == tid || face_tets[1] == tid);
    tid == face_tets[1]
}

impl TetSet {
    pub(crate) fn form_bbox(bbox: BBox, surfaces: &[Surf]) -> TetSet {
        const TETS: [[usize; 4]; 18] = [
            [0, 1, 7, 3],
            [7, 0, 5, 1],
            [4, 0, 5, 7],
            [4, 6, 0, 7],
            [0, 7, 6, 2],
            [7, 2, 0, 3],
            [1, 7, 3, INVALID_IND],
            [0, 1, 3, INVALID_IND],
            [1, 5, 7, INVALID_IND],
            [0, 5, 1, INVALID_IND],
            [0, 4, 5, INVALID_IND],
            [4, 7, 5, INVALID_IND],
            [0, 6, 4, INVALID_IND],
            [4, 6, 7, INVALID_IND],
            [0, 2, 6, INVALID_IND],
            [2, 7, 6, INVALID_IND],
            [0, 3, 2, INVALID_IND],
            [2, 3, 7, INVALID_IND],
        ];
        const TET_FACES: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];
        let points = vec![
            bbox.min[0],
            bbox.min[1],
            bbox.min[2],
            bbox.min[0],
            bbox.min[1],
            bbox.max[2],
            bbox.min[0],
            bbox.max[1],
            bbox.min[2],
            bbox.min[0],
            bbox.max[1],
            bbox.max[2],
            bbox.max[0],
            bbox.min[1],
            bbox.min[2],
            bbox.max[0],
            bbox.min[1],
            bbox.max[2],
            bbox.max[0],
            bbox.max[1],
            bbox.min[2],
            bbox.max[0],
            bbox.max[1],
            bbox.max[2],
        ];

        let hash_tri = |mut verts: [usize; 3]| {
            verts.sort();
            if verts[2] == INVALID_IND {
                verts[2] = 8;
            }
            (verts[2] << 6) | (verts[1] << 3) | verts[0]
        };

        let mut tet_faces = vec![[INVALID_IND; 4]; 18];
        let mut face_map = HashMap::<usize, usize>::with_capacity(36);

        let mut tet_face_vertices = Vec::with_capacity(36);
        for (tid, tet_vertices) in TETS.into_iter().enumerate() {
            for (i, face_indices) in TET_FACES.into_iter().enumerate() {
                let face_vertices = face_indices.map(|idx| tet_vertices[idx]);
                let fid = face_map.len();
                match face_map.entry(hash_tri(face_vertices)) {
                    hashbrown::hash_map::Entry::Occupied(entry) => {
                        tet_faces[tid][i] = twin_index(*entry.get());
                    }
                    hashbrown::hash_map::Entry::Vacant(entry) => {
                        let val = fid << 1;
                        entry.insert(val);
                        tet_faces[tid][i] = val;
                        tet_face_vertices.push(face_vertices);
                    }
                }
            }
        }

        let mut tet_mesh = SurfaceMesh::new(tet_face_vertices, std::alloc::Global);
        let mut face_tets = vec![[INVALID_IND; 2]; tet_mesh.n_faces()];
        let mut infinite_edges = [EdgeId::default(); 8];
        let square_edge_lengths = Vec::from_iter(tet_mesh.edges().map(|edge| {
            let eid = *edge;
            let [va, vb] = tet_mesh.e_vertices(eid);
            match [va.valid(), vb.valid()] {
                [_, false] => {
                    infinite_edges[va] = eid;
                    0.0
                }
                [false, _] => {
                    infinite_edges[vb] = eid;
                    0.0
                }
                _ => square_norm(&sub_short::<3, _>(
                    point::<3>(&points, *va),
                    point::<3>(&points, *vb),
                )),
            }
        }));
        let tets = tet_faces
            .into_iter()
            .enumerate()
            .zip(TETS)
            .map(|((tid, ori_faces), tet_vertices)| {
                let vertices = tet_vertices.map(VertexId::from);
                let mut edges = [EdgeId::default(); 6];
                for ((va, vb), edge) in vertices.into_iter().tuple_combinations().zip(&mut edges) {
                    *edge = match [va.valid(), vb.valid()] {
                        [_, false] => infinite_edges[va],
                        [false, _] => infinite_edges[vb],
                        _ => tet_mesh.e_from_va_vb(va, vb),
                    };
                }
                let faces = ori_faces.map(|ori_fid| {
                    let (fid, reversed) = decode_index(ori_fid);
                    if reversed {
                        face_tets[fid][1] = tid;
                    } else {
                        face_tets[fid][0] = tid;
                    }
                    fid.into()
                });
                let surface_evaluations = if tid < 6 {
                    let tet_points = tet_vertices.map(|idx| point::<3>(&points, idx));
                    TinyVec::from_iter(surfaces.iter().enumerate().map(|(sid, srf)| {
                        SurfaceEvaluation {
                            sid,
                            evaluation: tet_points.map(|p| srf.eval(p)),
                        }
                    }))
                } else {
                    TinyVec::new()
                };
                Tet {
                    vertices,
                    edges,
                    faces,
                    surface_evaluations,
                }
            })
            .collect_vec();

        let mesh_ptr = unsafe { std::mem::transmute::<_, *mut SurfaceMesh>(&mut tet_mesh) };
        for edge in tet_mesh.edges() {
            let mut tet_halfedges_map = HashMap::<usize, [HalfedgeId; 2]>::new();
            for he in edge.halfedges() {
                let hid = *he;
                let fid = *he.face();
                for tid in face_tets[fid] {
                    match tet_halfedges_map.entry(tid) {
                        hashbrown::hash_map::Entry::Occupied(mut entry) => {
                            entry.get_mut()[1] = hid;
                        }
                        hashbrown::hash_map::Entry::Vacant(entry) => {
                            entry.insert([hid, HalfedgeId::default()]);
                        }
                    }
                }
            }

            let first_hid = *edge.halfedge();
            let mut sorted_halfedges = vec![first_hid];
            let mut curr_hid = first_hid;
            let mut curr_tid = face_tets[tet_mesh.he_face(first_hid)][0];
            loop {
                curr_hid = {
                    let halfedges = tet_halfedges_map.get(&curr_tid).unwrap();
                    if halfedges[0] == curr_hid {
                        halfedges[1]
                    } else {
                        halfedges[0]
                    }
                };
                if curr_hid == first_hid {
                    break;
                }
                sorted_halfedges.push(curr_hid);

                let candidates = face_tets[tet_mesh.he_face(first_hid)];
                curr_tid = if candidates[0] == curr_tid {
                    candidates[1]
                } else {
                    candidates[0]
                };
            }
            for (h1, h2) in sorted_halfedges.into_iter().circular_tuple_windows() {
                unsafe { (*mesh_ptr).set_he_sibling(h1, h2) };
            }
        }

        TetSet {
            points,
            mesh: tet_mesh,
            tets,
            face_tets,
            square_edge_lengths,
        }
    }

    pub(crate) fn adaptive_subdivide(&mut self, srf_datum: &[SurfaceData], sq_eps: f64) {
        let mut pq = BinaryHeap::new();
        for tid in 0..6 {
            push_longest_edge(tid, self, srf_datum, &mut pq, sq_eps);
        }

        let split_bump = Bump::new();

        while !pq.is_empty() {
            let EdgeAndLen { eid, len } = pq.pop().unwrap();
            if self.square_edge_lengths[eid] != len {
                // edge changed
                continue;
            }
            let (new_vert, tet_pairs) = self.split_edge(eid, &split_bump);
            let p = point_3(&tets.points, new_vert.0);

            for (sid, surf) in data.surface_datum.iter().map(|d| &d.surf).enumerate() {
                data.vals_and_grads[sid].push(surf.eval(p));
            }

            for tid in tet_pairs.into_iter().flatten() {
                push_longest_edge(tid, tets, data, sq_eps);
            }
        }
    }

    /// return tet and its start face index
    pub(crate) fn tets_around_edge(&self, eid: EdgeId) -> TetsAroundEdge<'_> {
        TetsAroundEdge::new(self, eid)
    }

    fn edge_tets(&self, eid: EdgeId) -> impl Iterator<Item = (usize, FaceId)> {
        let first_hid = self.mesh.e_halfedge(eid);
    }

    pub(crate) fn build_descending_vertex_links(&self) -> Vec<VertexId> {
        let mut descent_links = vec![VertexId::default(); self.mesh.n_vertices()];
        for tet_verts in self.tets.iter().map(|tet| &tet.vertices) {
            let min_idx = tet_verts
                .map(|vid| point_3(&self.points, vid.0))
                .into_iter()
                .enumerate()
                .min_by(|&(_, pa), &(_, pb)| pa.partial_cmp(&pb).unwrap())
                .unwrap()
                .0;
            let min_v = tet_verts[min_idx];

            for &vid in tet_verts {
                if vid == min_v {
                    continue;
                }
                let next_vid = descent_links[vid];
                if !next_vid.valid() {
                    descent_links[vid] = min_v;
                }
            }
        }
        descent_links
    }

    pub(crate) fn split_edge<A: Allocator + Copy>(
        &mut self,
        split_eid: EdgeId,
        alloc: A,
    ) -> (VertexId, Vec<[usize; 2], A>) {
        let [va, vb] = self.mesh.e_vertices(split_eid);
        let new_vid = self.mesh.split_edge(split_eid, alloc);
        let new_eid = self.mesh.he_prev(self.mesh.e_halfedge(split_eid));
        {
            // update points
            let pa = point_3(&self.points, va.0);
            let pb = point_3(&self.points, vb.0);
            self.points.extend_from_slice(&[
                (pa[0] + pb[0]) * 0.5,
                (pa[1] + pb[1]) * 0.5,
                (pa[2] + pb[2]) * 0.5,
            ]);
        }
        {
            // update edge lengths
            let new_edge_square_len = self.square_edge_lengths[split_eid] * 0.25;
            self.square_edge_lengths.push(new_edge_square_len);
            self.square_edge_lengths[split_eid] = new_edge_square_len;
        }

        let mut side_halfedges = Vec::with_capacity_in(4, alloc);
        unsafe {
            let mesh_ptr = std::mem::transmute::<_, *mut SurfaceMesh>(&mut self.mesh);
            side_halfedges.extend(self.mesh.edge(split_eid).halfedges().map(|he| {
                let vc = *he.next().to();
                let fid = *he.face();
                let hid = (*mesh_ptr).split_face(fid, new_vid, vc);
                if self.mesh.he_to(hid) == vc {
                    [hid, self.mesh.he_twin(hid)]
                } else {
                    [self.mesh.he_sibling(hid), hid]
                }
            }));
        };

        let mut prev_tid = INVALID_IND;
        for (left_halfedges, right_halfedges) in side_halfedges.iter().circular_tuple_windows() {
            let [fa, fb] = [left_halfedges, right_halfedges].map(|halfedges| {
                let [f1, f2] = [self.mesh.he_face(halfedges[0]), self.mesh.he_face(halfedges[1])];
                if *f1 < *f2 {
                    f1
                } else {
                    f2
                }
            });

            let old_tid = if prev_tid == INVALID_IND {
                (|| {
                    for tid in self.face_tets[fa] {
                        if self.face_tets[fb].contains(&tid) {
                        return tid;
                        }
                    }
                    INVALID_IND
                })()
            } else {
                let face_tets = self.face_tets[fa];
                if face_tets[0] == prev_tid {
                    face_tets[1]
                } else {
                    face_tets[0]
                }
            };
            debug_assert!(old_tid != INVALID_IND);
            prev_tid = old_tid;

            let tet = &self.tets[old_tid];
            let [bottom_fid, top_fid] = tet.face_from_edge(split_eid, va);
            let bottom_hid = self.mesh.face(bottom_fid).halfedges().find_map(|he| {
                if *he.next().to() == va {
                    Some(*he)
                } else {
                    None
                }
            }).unwrap();
            let ori_bottom_hid = if self.mesh.he_to(bottom_hid) != self.mesh.he_to(right_halfedges[0]) {
                self.mesh.he_twin(bottom_hid)
            } else {
                bottom_hid
            };

            debug_assert!(self.mesh.he_to(left_halfedges[0]) == self.mesh.he_from(bottom_hid));
            debug_assert!(self.mesh.he_to(bottom_hid) == self.mesh.he_to(right_halfedges[0]));

            let new_fid = self.mesh.add_face_by_halfedges(&[left_halfedges[0], ori_bottom_hid, right_halfedges[1]], false);
            self.face_tets.push([self.tets.len(), old_tid]);
            // set sibling halfedges
            {
                let mut new_hid = self.mesh.f_halfedge(new_fid);
                debug_assert!(self.mesh.he_edge(new_hid) == self.mesh.he_edge(left_halfedges[0]));
                if self.mesh.he_sibling(left_halfedges[0]) == left_halfedges[1] {
                    self.mesh.insert_he_sibling(left_halfedges[0], new_hid);
                } else {
                    debug_assert!(self.mesh.he_sibling(left_halfedges[1]) == left_halfedges[0]);
                    self.mesh.insert_he_sibling(left_halfedges[1], new_hid);
                }

                new_hid = self.mesh.he_next(new_hid);
                debug_assert!(self.mesh.he_edge(new_hid) == self.mesh.he_edge(bottom_hid));
                let oppo_bottom_hid = self.mesh.face(top_fid).halfedges().find_map(|he| {
                    if *he.next().to() == vb {
                        Some(*he)
                    } else {
                        None
                    }
                }).unwrap();
                if self.mesh.he_sibling(bottom_hid) == oppo_bottom_hid {
                    self.mesh.insert_he_sibling(bottom_hid, new_hid);
                } else {
                    debug_assert!(self.mesh.he_sibling(oppo_bottom_hid) == bottom_hid);
                    self.mesh.insert_he_sibling(oppo_bottom_hid, new_hid);
                }

                new_hid = self.mesh.he_next(new_hid);
                debug_assert!(self.mesh.he_edge(new_hid) == self.mesh.he_edge(right_halfedges[0]));
                if self.mesh.he_sibling(right_halfedges[0]) == right_halfedges[1] {
                    self.mesh.insert_he_sibling(right_halfedges[0], new_hid);
                } else {
                    debug_assert!(self.mesh.he_sibling(right_halfedges[1]) == right_halfedges[0]);
                    self.mesh.insert_he_sibling(right_halfedges[1], new_hid);
                }
            }

            let [vc, vd] = [self.mesh.he_to(left_halfedges[0]), self.mesh.he_to(right_halfedges[0])];
            let [left_bottom_fid, left_top_fid] = if *self.mesh.halfedge(left_halfedges[0]).next().to() == va {
                debug_assert!(*self.mesh.halfedge(left_halfedges[1]).next().to() == vb);
                [self.mesh.he_face(left_halfedges[0]), self.mesh.he_face(left_halfedges[1])]
            } else {
                debug_assert!(*self.mesh.halfedge(left_halfedges[0]).next().to() == vb);
                debug_assert!(*self.mesh.halfedge(left_halfedges[1]).next().to() == va);
                [self.mesh.he_face(left_halfedges[1]), self.mesh.he_face(left_halfedges[0])]
            };

            let [right_bottom_fid, right_top_fid] = if *self.mesh.halfedge(right_halfedges[0]).next().to() == va {
                debug_assert!(*self.mesh.halfedge(right_halfedges[1]).next().to() == vb);
                [self.mesh.he_face(right_halfedges[0]), self.mesh.he_face(right_halfedges[1])]
            } else {
                debug_assert!(*self.mesh.halfedge(right_halfedges[0]).next().to() == vb);
                debug_assert!(*self.mesh.halfedge(right_halfedges[1]).next().to() == va);
                [self.mesh.he_face(right_halfedges[1]), self.mesh.he_face(right_halfedges[0])]
            };

            if !vc.valid() {
                let new_tet = {
                    Tet {
                        vertices: [vb, new_vid, vd, vc],
                        edges: [new_eid, self.mesh.he_edge(right_halfedges[0]), ]
                    }
                };
            }
        }

        /*let [va, vb] = self.mesh.e_vertices(eid);
        let mut tet_faces_map: hashbrown::HashMap<usize, [usize; 2], _, _> =
            hashbrown::HashMap::<usize, [usize; 2], _, _>::new_in(alloc);
        let mut faces = Vec::new_in(alloc);
        let mut oppo_verts = Vec::new_in(alloc);
        for (i, he) in self.mesh.edge(eid).halfedges().enumerate() {
            let fid = self.mesh.he_face(*he);
            faces.push(fid);
            for tid in self.face_tets[fid.0].iter() {
                if tid == &INVALID_IND {
                    continue;
                }
                if let Some(tet_faces) = tet_faces_map.get_mut(tid) {
                    tet_faces[1] = i;
                } else {
                    tet_faces_map.insert(*tid, [i, INVALID_IND]);
                }
            }
            oppo_verts.push(*he.next().to());
        }
        let mut tet_faces = Vec::from_iter(tet_faces_map.clone());
        tet_faces.sort_unstable();
        let mut oppo_halfedges = Vec::with_capacity_in(tet_faces_map.len(), alloc);
        let mut bottom_faces = Vec::with_capacity_in(tet_faces_map.len(), alloc);
        let mut top_faces = Vec::with_capacity_in(tet_faces_map.len(), alloc);
        for &(tid, [fa, fb]) in &tet_faces {
            let bottom_top_faces = || {
                let mut fc = FaceId::default();
                let mut fd = FaceId::default();
                for fid in self.tets[tid].faces {
                    if fid == faces[fa] || fid == faces[fb] {
                        continue;
                    }
                    if fc.0 == INVALID_IND {
                        fc = fid;
                    } else {
                        fd = fid;
                    }
                }
                let h = self.mesh.fv_halfedge(fc, va);
                if h.0 != INVALID_IND {
                    [fc, fd]
                } else {
                    [fd, fc]
                }
            };
            let [bottom_fid, top_fid] = bottom_top_faces();
            debug_assert!(self.mesh.fv_halfedge(bottom_fid, va).0 != INVALID_IND);
            debug_assert!(self.mesh.fv_halfedge(top_fid, vb).0 != INVALID_IND);

            bottom_faces.push(bottom_fid);
            top_faces.push(top_fid);
            let hid = self.mesh.fv_halfedge(bottom_fid, va);

            if tet_face_reversed(&self.face_tets[bottom_fid], tid) {
                oppo_halfedges.push(self.mesh.he_twin(self.mesh.he_next(hid)));
            } else {
                oppo_halfedges.push(self.mesh.he_next(hid));
            }
        }
        let new_vert = self.mesh.split_edge(eid, alloc);
        let [bottom_eid, top_eid] = {
            let hid = self.mesh.e_halfedge(eid);
            if self.mesh.he_from(hid) == va || self.mesh.he_to(hid) == va {
                [eid, self.square_edge_lengths.len().into()]
            } else {
                [self.square_edge_lengths.len().into(), eid]
            }
        };
        {
            // update points
            let pa = point_3(&self.points, va.0);
            let pb = point_3(&self.points, vb.0);
            self.points.extend_from_slice(&[
                (pa[0] + pb[0]) * 0.5,
                (pa[1] + pb[1]) * 0.5,
                (pa[2] + pb[2]) * 0.5,
            ]);
        }
        {
            // update edge lengths
            let new_edge_square_len = self.square_edge_lengths[eid] * 0.25;
            self.square_edge_lengths.push(new_edge_square_len);
            self.square_edge_lengths[eid] = new_edge_square_len;
        }

        let faces_capacity = self.mesh.n_faces() + faces.len() + tet_faces_map.len();
        self.face_tets.reserve(faces_capacity);
        let mut new_halfedges = Vec::with_capacity_in(oppo_verts.len(), alloc);
        new_halfedges.extend(oppo_verts.into_iter().zip(&faces).map(|(v, &fid)| {
            let hid = self.mesh.split_face(fid, new_vert, v);
            self.square_edge_lengths.push(square_edge_length(
                &self.points,
                self.mesh.he_edge(hid),
                &self.mesh,
            ));
            self.face_tets.push(self.face_tets[fid]);
            if self.mesh.he_to(hid) != new_vert {
                hid
            } else {
                self.mesh.he_twin(hid)
            }
        }));

        debug_assert!(
            new_halfedges
                .iter()
                .all(|&hid| self.mesh.he_from(hid) == new_vert)
        );

        let mut result_tets = Vec::with_capacity_in(tet_faces_map.len() << 1, alloc);
        for ((((tid, face_indices), oppo_hid), bottom_fid), top_fid) in tet_faces
            .into_iter()
            .zip(oppo_halfedges)
            .zip(bottom_faces)
            .zip(top_faces)
        {
            let mut ha = new_halfedges[face_indices[0]];
            let mut hb = new_halfedges[face_indices[1]];
            if self.mesh.he_to(ha) != self.mesh.he_from(oppo_hid) {
                std::mem::swap(&mut ha, &mut hb);
            }
            hb = self.mesh.he_twin(hb);

            debug_assert!(self.mesh.he_to(ha) == self.mesh.he_from(oppo_hid));
            debug_assert!(self.mesh.he_from(hb) == self.mesh.he_to(oppo_hid));
            debug_assert!(self.mesh.he_to(hb) == self.mesh.he_from(ha));

            let new_fid = self.mesh.add_face_by_halfedges(&[ha, oppo_hid, hb]);
            let new_tid = self.tets.len();

            // We have chose `oppo_hid` corresponding to bottom face with positive orientation
            self.face_tets.push([tid, new_tid]);

            let replace = |tets: &mut [usize; 2]| {
                for t in tets {
                    if *t == tid {
                        *t = new_tid;
                        break;
                    }
                }
            };

            replace(&mut self.face_tets[bottom_fid]);

            let side_faces_and_edges = [ha, self.mesh.he_twin(hb)].map(|hid| {
                // when the face is split, the f_halfedge of the face is the split halfedge
                let he = self.mesh.halfedge(hid);
                let next_he = he.next();
                if *next_he.to() == va {
                    let fid = self.mesh.he_face(*he);
                    replace(&mut self.face_tets[fid]);

                    let prev_twin = he.twin().prev();

                    debug_assert!(*prev_twin.from() == vb);
                    debug_assert!(*prev_twin.to() == *he.to());

                    (
                        [self.mesh.he_face(*prev_twin), fid],
                        [*prev_twin.edge(), *next_he.edge()],
                    )
                } else {
                    debug_assert!(*next_he.to() == vb);
                    let prev_twin = he.twin().prev();
                    debug_assert!(*prev_twin.from() == va);
                    debug_assert!(*prev_twin.to() == *he.to());
                    let twin_fid = self.mesh.he_face(*prev_twin);
                    replace(&mut self.face_tets[twin_fid]);
                    (
                        [self.mesh.he_face(*he), twin_fid],
                        [*next_he.edge(), *prev_twin.edge()],
                    )
                }
            });

            let [vc, vd] = self.mesh.he_vertices(oppo_hid);
            let [ea, eb, oppo_eid] = [ha, hb, oppo_hid].map(|hid| self.mesh.he_edge(hid));

            self.tet_vertices[tid] = [new_vert, vb, vd, vc];
            self.tet_edges[tid] = [
                top_eid,
                eb,
                ea,
                side_faces_and_edges[1].1[0],
                side_faces_and_edges[0].1[0],
                oppo_eid,
            ];
            self.tet_faces[tid] = [
                top_fid,
                new_fid,
                side_faces_and_edges[0].0[0],
                side_faces_and_edges[1].0[0],
            ];

            self.tet_vertices.push([new_vert, va, vc, vd]);
            self.tet_edges.push([
                bottom_eid,
                ea,
                eb,
                side_faces_and_edges[0].1[1],
                side_faces_and_edges[1].1[1],
                oppo_eid,
            ]);

            self.tet_faces.push([
                bottom_fid,
                new_fid,
                side_faces_and_edges[1].0[1],
                side_faces_and_edges[0].0[1],
            ]);
            self.surf_indices.push(self.surf_indices[tid].clone());
            #[cfg(debug_assertions)]
            {
                for _t in [tid, new_tid] {
                    let _verts = &self.tet_vertices[_t];
                    let _pts = _verts.map(|vid| point_3(&self.points, vid.0));
                    let _ori = crate::predicates::orient3d::orient3d_eeee(
                        &_pts[0], &_pts[1], &_pts[2], &_pts[3], alloc,
                    );
                    debug_assert!(_ori.is_pos());
                }
            }

            result_tets.push([tid, new_tid]);
        }
        (new_vert, result_tets)*/
        (VertexId::default(), Vec::new_in(alloc))
    }
}

fn push_longest_edge(
    tid: usize,
    tets: &mut TetSet,
    srf_datum: &[SurfaceData],
    pq: &mut BinaryHeap<EdgeAndLen>,
    sq_eps: f64,
) {
    let tet = &mut tets.tets[tid];
    if tet.subdividable(&tets.points, &tets.square_edge_lengths, &srf_datum, sq_eps) {
        let longest_eid = *tet
            .edges
            .iter()
            .max_by(|&&ea, &&eb| {
                tets.square_edge_lengths[ea]
                    .partial_cmp(&tets.square_edge_lengths[eb])
                    .unwrap()
            })
            .unwrap();
        pq.push(EdgeAndLen {
            eid: longest_eid,
            len: tets.square_edge_lengths[longest_eid],
        });
    }
}

pub(crate) struct TetsAroundEdge<'a> {
    eid: EdgeId,
    tid: usize,
    face_index_in_tet: usize,
    first_tid: usize,
    curr_fid: FaceId,
    next_fid: FaceId,
    is_first: bool,
    tets: &'a TetSet,
}

pub(crate) const EDGE_FACE_INDICES: [[usize; 2]; 6] =
    [[2, 3], [1, 3], [1, 2], [0, 3], [0, 2], [0, 1]];

impl<'a> TetsAroundEdge<'a> {
    fn new(tets: &'a TetSet, eid: EdgeId) -> Self {
        let mesh = &tets.mesh;
        let curr_fid = mesh
            .edge(eid)
            .halfedges()
            .map(|he| *he.face())
            .find(|&fid| tets.face_tets[fid][1] == INVALID_IND)
            .unwrap_or(*mesh.edge(eid).halfedge().face());
        let tid = tets.face_tets[curr_fid][0];

        let mut data = Self {
            eid,
            tid,
            face_index_in_tet: 0,
            first_tid: tid,
            curr_fid,
            next_fid: FaceId::default(),
            is_first: true,
            tets,
        };
        data.compute_next_face();
        data
    }

    fn compute_next_face(&mut self) {
        // let tets = &self.tets;
        // let edge_index = tets.tet_edges[self.tid]
        //     .iter()
        //     .position(|&eid| self.eid == eid)
        //     .unwrap();

        // let curr_face_indices = EDGE_FACE_INDICES[edge_index];
        // let tet_faces = &tets.tet_faces[self.tid];
        // let curr_face_pos = curr_face_indices
        //     .into_iter()
        //     .position(|idx| tet_faces[idx] == self.curr_fid)
        //     .unwrap();

        // self.face_index_in_tet = curr_face_indices[curr_face_pos];
        // self.next_fid = tet_faces[curr_face_indices[curr_face_pos ^ 1]];
    }

    fn is_end(&self) -> bool {
        self.tid == INVALID_IND || (!self.is_first && self.tid == self.first_tid)
    }

    fn next_self(&mut self) {
        self.curr_fid = self.next_fid;

        let tets = &self.tets;
        let cells = &tets.face_tets[self.curr_fid];
        debug_assert!(cells[0] == self.tid || cells[1] == self.tid);

        self.tid = if cells[0] == self.tid {
            cells[1]
        } else {
            cells[0]
        };

        if self.is_end() {
            return;
        }

        self.compute_next_face();
    }
}

impl<'a> Iterator for TetsAroundEdge<'a> {
    type Item = [usize; 2];

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_end() {
            return None;
        }
        let ret = Some([self.tid, self.face_index_in_tet]);
        self.next_self();
        ret
    }
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
        let p = point::<2>(&points, vid);
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
                let p = point::<3>(&points, vid);
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
