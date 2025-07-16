use std::{
    alloc::Allocator,
    collections::BinaryHeap,
};

use bumpalo::Bump;
use hashbrown::{HashMap, HashSet};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    INVALID_IND, decode_index,
    geometry::{BBox, Surf, Surface},
    math::{cross, cross_in, dot, square_norm, sub_short},
    mesh::{EdgeId, ElementId, FaceId, HalfedgeId, Mesh, SurfaceMesh, VertexId},
    point, point_3,
    triangle::{convex_2, convex_3},
    twin_index,
};

#[derive(Default, Clone)]
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
    pub(crate) fn valid(&self) -> bool {
        self.vertices[3].valid()
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
        if !self.valid()
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
            let mut sorted_halfedges = Vec::with_capacity(tet_halfedges_map.len());
            sorted_halfedges.push(first_hid);
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

                let candidates = face_tets[tet_mesh.he_face(curr_hid)];
                curr_tid = if candidates[0] == curr_tid {
                    candidates[1]
                } else {
                    candidates[0]
                };
            }

            debug_assert!(sorted_halfedges.len() == tet_halfedges_map.len());
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

        let mut split_bump = Bump::new();

        while !pq.is_empty() {
            let EdgeAndLen { eid, len } = pq.pop().unwrap();
            if self.square_edge_lengths[eid] != len {
                // edge changed
                continue;
            }
            split_bump.reset();

            let start_face_idx =
                self.split_edge(eid, srf_datum, &split_bump);

            for fid in start_face_idx..self.face_tets.len() {
                for tid in self.face_tets[fid] {
                    push_longest_edge(tid, self, srf_datum, &mut pq, sq_eps);
                }
            }
        }
    }

    /// return tet and its start face index
    pub(crate) fn tets_around_edge(&self, eid: EdgeId) -> TetsAroundEdge<'_> {
        TetsAroundEdge::new(self, eid)
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
        srf_datum: &[SurfaceData],
        alloc: A,
    ) -> usize {
        let [va, vb] = self.mesh.e_vertices(split_eid);
        let new_eid = self.mesh.n_edges_capacity().into();
        let ve = self.mesh.split_edge(split_eid, alloc);
        let new_pt = {
            // update points
            let pa = point_3(&self.points, va.0);
            let pb = point_3(&self.points, vb.0);
            self.points.extend_from_slice(&[
                (pa[0] + pb[0]) * 0.5,
                (pa[1] + pb[1]) * 0.5,
                (pa[2] + pb[2]) * 0.5,
            ]);
            point::<3>(&self.points, ve.0)
        };
        {
            // update edge lengths
            let new_edge_square_len = self.square_edge_lengths[split_eid] * 0.25;
            self.square_edge_lengths.push(new_edge_square_len);
            self.square_edge_lengths[split_eid] = new_edge_square_len;
        }
        let mut side_halfedges =
            Vec::with_capacity_in(self.mesh.edge(split_eid).halfedges().count(), alloc);
        let mesh_ptr = unsafe { std::mem::transmute::<_, *mut SurfaceMesh>(&mut self.mesh) };
        let mesh_ref = &self.mesh; // Cache mesh reference to reduce dereferencing
        side_halfedges.extend(mesh_ref.edge(split_eid).halfedges().map(|he| {
            let ([h_ac, h_bc], vc) = if *he.to() == vb {
                let he_next = he.next();
                ([*he.prev().prev(), *he_next], *he_next.to())
            } else {
                let he_prev = he.prev();
                ([*he.next().next(), *he.prev()], *he_prev.from())
            };

            let fid = *he.face();
            let new_hid = unsafe { (*mesh_ptr).split_face(fid, ve, vc) };
            self.face_tets.push(self.face_tets[fid]);
            if vc.valid() {
                self.square_edge_lengths
                    .push(square_norm(&sub_short::<3, _>(
                        point::<3>(&self.points, *vc),
                        new_pt,
                    )));
            } else {
                self.square_edge_lengths.push(0.0);
            }

            debug_assert!(mesh_ref.he_to(new_hid) == vc);
            debug_assert!(mesh_ref.he_from(new_hid) == ve);
            [new_hid, mesh_ref.he_sibling(new_hid), h_ac, h_bc]
        }));

        let ret = self.face_tets.len();
        let mut prev_tid = INVALID_IND;
        let mut evaluation_map = HashMap::with_capacity_in(
            self.tets[self.face_tets[*mesh_ref.edge(split_eid).halfedge().face()][0]]
                .surface_evaluations
                .len(),
            alloc,
        );

        for ([h_ec, h_ce, h_ac, h_bc], [h_ed, h_de, h_ad, h_bd]) in
            side_halfedges.into_iter().circular_tuple_windows()
        {
            let [
                left_bottom_fid,
                left_top_fid,
                right_bottom_fid,
                right_top_fid,
            ] = [
                *mesh_ref.halfedge(h_ac).face(),
                *mesh_ref.halfedge(h_bc).face(),
                *mesh_ref.halfedge(h_ad).face(),
                *mesh_ref.halfedge(h_bd).face(),
            ];
            let [vc, vd] = [mesh_ref.he_to(h_ec), mesh_ref.he_to(h_ed)];

            // Optimized tetrahedron finding
            let old_tid = if prev_tid == INVALID_IND {
                let left_tets = self.face_tets[left_bottom_fid];
                let right_tets = self.face_tets[right_bottom_fid];
                if left_tets[0] == right_tets[0] || left_tets[0] == right_tets[1] {
                    left_tets[0]
                } else {
                    left_tets[1]
                }
            } else {
                let face_tets = self.face_tets[left_bottom_fid];
                if face_tets[0] == prev_tid {
                    face_tets[1]
                } else {
                    face_tets[0]
                }
            };
            debug_assert!(old_tid != INVALID_IND);
            prev_tid = old_tid;
            let new_tid = self.tets.len();
            let tet = &mut self.tets[old_tid];

            // Get immutable reference first to avoid borrowing conflicts
            let [bottom_fid, top_fid] = tet.face_from_edge(split_eid, va);

            let bottom_hid = mesh_ref.get_he_from_oppo_vertex(bottom_fid, va);

            let h_cd = if mesh_ref.he_to(bottom_hid) != vd {
                mesh_ref.he_twin(bottom_hid)
            } else {
                bottom_hid
            };

            debug_assert!(mesh_ref.he_from(h_cd) == mesh_ref.he_to(h_ec));
            debug_assert!(mesh_ref.he_to(h_cd) == mesh_ref.he_to(h_ed));

            let new_fid = unsafe { (*mesh_ptr).add_face_by_halfedges(&[h_ec, h_cd, h_de], false) };
            self.face_tets.push([new_tid, old_tid]);

            // set sibling halfedges - optimized by caching halfedge lookups
            {
                let set_sibling = |prev_hid, next_hid, curr_hid| {
                    debug_assert!(mesh_ref.he_edge(curr_hid) == mesh_ref.he_edge(prev_hid));
                    debug_assert!(mesh_ref.he_edge(curr_hid) == mesh_ref.he_edge(next_hid));
                    if mesh_ref.he_sibling(prev_hid) == next_hid {
                        unsafe {
                            (*mesh_ptr).insert_he_sibling(prev_hid, curr_hid);
                        }
                    } else {
                        debug_assert!(mesh_ref.he_sibling(next_hid) == prev_hid);
                        unsafe {
                            (*mesh_ptr).insert_he_sibling(next_hid, curr_hid);
                        }
                    }
                };

                let mut new_hid = mesh_ref.f_halfedge(new_fid);
                set_sibling(h_ec, h_ce, new_hid);

                new_hid = mesh_ref.he_next(new_hid);
                // Optimized opposite halfedge finding
                let oppo_bottom_hid = mesh_ref
                    .face(top_fid)
                    .halfedges()
                    .find_map(|he| {
                        if *he.next().to() == vb {
                            Some(*he)
                        } else {
                            None
                        }
                    })
                    .unwrap();
                set_sibling(bottom_hid, oppo_bottom_hid, new_hid);

                new_hid = mesh_ref.he_next(new_hid);
                set_sibling(h_ed, h_de, new_hid);
            }

            // Optimized vertex index mapping using lookup table instead of loop
            let mut tet_indices = [usize::MAX; 4]; // Use MAX as sentinel
            for (i, &v) in tet.vertices.iter().enumerate() {
                match v {
                    v if v == va => tet_indices[0] = i,
                    v if v == vb => tet_indices[1] = i,
                    v if v == vc => tet_indices[2] = i,
                    v if v == vd => tet_indices[3] = i,
                    _ => {} // Should not happen for valid tetrahedron
                }
            }
            let mut get_srf_eval = |indices: [usize; 4]| {
                TinyVec::from_iter(tet.surface_evaluations.iter().map(|tet_eval| {
                    let sid = tet_eval.sid;
                    let evaluation = indices.map(|idx| {
                        if idx < 4 {
                            tet_eval.evaluation[tet_indices[idx]]
                        } else {
                            *evaluation_map
                                .entry(sid)
                                .or_insert_with(|| srf_datum[sid].surf.eval(new_pt))
                        }
                    });
                    SurfaceEvaluation { sid, evaluation }
                }))
            };

            // Batch edge lookups to reduce function calls
            let [e_ec, e_ed, e_cd, e_ac, e_bc, e_ad, e_bd] = [
                mesh_ref.he_edge(h_ec), // e_ec
                mesh_ref.he_edge(h_ed), // e_ed
                mesh_ref.he_edge(h_cd), // e_cd
                mesh_ref.he_edge(h_ac), // e_ac
                mesh_ref.he_edge(h_bc), // e_bc
                mesh_ref.he_edge(h_ad), // e_ad
                mesh_ref.he_edge(h_bd), // e_bd
            ];

            // Batch face updates with early termination
            for fid in [left_top_fid, right_top_fid, top_fid] {
                let face_tet_ref = &mut self.face_tets[fid];
                if face_tet_ref[0] == old_tid {
                    face_tet_ref[0] = new_tid;
                } else {
                    debug_assert!(face_tet_ref[1] == old_tid);
                    face_tet_ref[1] = new_tid;
                }
            }

            let new_tet = if !vc.valid() {
                let new_tet = Tet {
                    vertices: [vb, ve, vd, vc],
                    edges: [split_eid, e_bd, e_bc, e_ed, e_ec, e_cd],
                    faces: [new_fid, top_fid, left_top_fid, right_top_fid],
                    surface_evaluations: get_srf_eval([1, 4, 3, 2]),
                };
                tet.vertices = [ve, va, vd, vc];
                tet.edges = [new_eid, e_ed, e_ec, e_ad, e_ac, e_cd];
                tet.faces = [bottom_fid, new_fid, left_bottom_fid, right_bottom_fid];
                tet.surface_evaluations = get_srf_eval([4, 0, 3, 2]);
                new_tet
            } else {
                let new_tet = Tet {
                    vertices: [ve, vb, vc, vd],
                    edges: [split_eid, e_ec, e_ed, e_bc, e_bd, e_cd],
                    faces: [top_fid, new_fid, right_top_fid, left_top_fid],
                    surface_evaluations: get_srf_eval([4, 1, 2, 3]),
                };
                tet.vertices = [va, ve, vc, vd];
                tet.edges = [new_eid, e_ac, e_ad, e_ec, e_ed, e_cd];
                tet.faces = [new_fid, bottom_fid, right_bottom_fid, left_bottom_fid];
                tet.surface_evaluations = get_srf_eval([0, 4, 2, 3]);
                new_tet
            };

            self.tets.push(new_tet);

            #[cfg(debug_assertions)]
            {
                let check_tet = |tid: usize| {
                    let tet = &self.tets[tid];
                    for ((va, vb), eid) in
                        tet.vertices.into_iter().tuple_combinations().zip(tet.edges)
                    {
                        let [vc, vd] = mesh_ref.e_vertices(eid);
                        debug_assert!((va == vc && vb == vd) || (va == vd && vb == vc));
                        if va.valid() && vb.valid() {
                            let pa = point::<3>(&self.points, *va);
                            let pb = point::<3>(&self.points, *vb);
                            let sq_len = square_norm(&sub_short::<3, _>(pa, pb));
                            debug_assert!((sq_len - self.square_edge_lengths[eid]).abs() < 1e-12);
                        }

                        let halfedges = Vec::from_iter(mesh_ref.edge(eid).halfedges());
                        for (he1, he2) in halfedges.iter().circular_tuple_windows() {
                            let f1 = *he1.face();
                            let f2 = *he2.face();
                            let [t1, t2] = self.face_tets[f1];
                            let [t3, t4] = self.face_tets[f2];
                            assert!(t1 == t3 || t1 == t4 || t2 == t3 || t2 == t4);
                        }
                    }

                    if tet.vertices[3].valid() {
                        let points = tet.vertices.map(|vid| point::<3>(&self.points, *vid));
                        for eval in &tet.surface_evaluations {
                            for i in 0..4 {
                                let e1 = eval.evaluation[i];
                                let e2 = srf_datum[eval.sid].surf.eval(points[i]);
                                let length = square_norm(&sub_short::<4, _>(&e1, &e2));
                                debug_assert!(length < 1e-12);
                            }
                        }
                    }

                    let mut tet_vertices = tet.vertices.map(|vid| *vid);
                    tet_vertices.sort();
                    for (fid, vd) in tet.faces.into_iter().zip(tet.vertices) {
                        let he = mesh_ref.face(fid).halfedge();
                        let mut vs0 = [**he.from(), **he.to(), **he.next().to(), *vd];
                        vs0.sort();
                        debug_assert!(vs0 == tet_vertices);
                    }
                };

                check_tet(old_tid);
                check_tet(new_tid);
            }
        }

        ret
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
