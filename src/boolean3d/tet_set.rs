use std::{alloc::Allocator, collections::BinaryHeap, ptr::NonNull};

use bumpalo::Bump;
use hashbrown::{HashMap, HashSet};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    INVALID_IND, decode_index,
    geometry::{BBox, Surf, Surface},
    is_positive,
    math::{cross, cross_in, dot, square_norm, sub_short},
    mesh1::{
        EdgeHalfedge, EdgeId, ElementId, FaceId, FaceMut, HalfedgeDataExt, HalfedgeId,
        HalfedgeNavigation, HalfedgeNavigationMut, Mesh, MeshCore, SurfaceMesh, VertexId,
    },
    point,
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
}

#[derive(Default, Clone)]
pub(crate) struct VertexData {
    pt: [f64; 3],
}

#[derive(Default, Clone)]
pub(crate) struct EdgeData {
    square_len: f64,
}

#[derive(Default, Clone)]
pub(crate) struct FaceData {
    tets: [usize; 2],
}

pub(crate) struct TetSet {
    pub(crate) mesh: SurfaceMesh<VertexData, (), EdgeData, FaceData, std::alloc::Global>,
    pub(crate) tets: Vec<Tet>,
}

#[inline]
pub(crate) fn tet_face_reversed(face_tets: &[usize; 2], tid: usize) -> bool {
    debug_assert!(face_tets[0] == tid || face_tets[1] == tid);
    tid == face_tets[1]
}

impl TetSet {
    pub(crate) fn from_bbox(bbox: BBox, surfaces: &[Surf]) -> TetSet {
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

        let mut tet_mesh =
            SurfaceMesh::<VertexData, (), EdgeData, FaceData, std::alloc::Global>::new_in(
                tet_face_vertices,
                std::alloc::Global,
            );
        let mut mesh_ptr = NonNull::from_ref(&tet_mesh);
        tet_mesh.vertices_mut().for_each(|vertex| {
            let id = *vertex.id;
            let pt = &mut vertex.data.property.pt;
            pt[0] = if id < 4 { bbox.min[0] } else { bbox.max[0] };
            pt[1] = if ((id % 4) >> 1) == 0 {
                bbox.min[1]
            } else {
                bbox.max[1]
            };
            pt[2] = if is_positive(id) {
                bbox.min[2]
            } else {
                bbox.max[2]
            };
        });

        let all_surfaces_evaluations = surfaces
            .iter()
            .map(|srf| {
                Vec::from_iter(
                    tet_mesh
                        .vertex_datum()
                        .map(|data| srf.eval(&data.property.pt)),
                )
            })
            .collect::<Vec<_>>();

        let mut infinite_edges = [EdgeId::default(); 8];
        unsafe {
            mesh_ptr.as_mut().edges_mut().for_each(|edge| {
                let eid = edge.id;
                let [va, vb] = tet_mesh.e_vertices(eid);
                match [va.valid(), vb.valid()] {
                    [_, false] => {
                        infinite_edges[*va] = eid;
                    }
                    [false, _] => {
                        infinite_edges[*vb] = eid;
                    }
                    _ => {
                        let pa = &tet_mesh.vertex(va).data.property.pt;
                        let pb = &tet_mesh.vertex(vb).data.property.pt;
                        edge.data.property.square_len = square_norm(&sub_short::<3, _>(pa, pb));
                    }
                }
            });
        }
        let tets = tet_faces
            .into_iter()
            .enumerate()
            .zip(TETS)
            .map(|((tid, ori_faces), tet_vertices)| {
                let vertices = tet_vertices.map(VertexId::from);
                let mut edges = [EdgeId::default(); 6];
                for ((va, vb), edge) in vertices.into_iter().tuple_combinations().zip(&mut edges) {
                    *edge = match [va.valid(), vb.valid()] {
                        [_, false] => infinite_edges[*va],
                        [false, _] => infinite_edges[*vb],
                        _ => tet_mesh.e_from_vertices(va, vb),
                    };
                }
                let faces = ori_faces.map(|ori_fid| {
                    let (fid, reversed) = decode_index(ori_fid);
                    if reversed {
                        tet_mesh.face_data_mut(fid.into()).property.tets[1] = tid;
                    } else {
                        tet_mesh.face_data_mut(fid.into()).property.tets[0] = tid;
                    }
                    fid.into()
                });
                let surface_evaluations = if tid < 6 {
                    TinyVec::from_iter((0..surfaces.len()).map(|sid| SurfaceEvaluation {
                        sid,
                        evaluation: tet_vertices.map(|idx| all_surfaces_evaluations[sid][idx]),
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

        for edge in tet_mesh.edges() {
            let mut tet_halfedges_map = HashMap::<usize, [HalfedgeId; 2]>::new();
            for he in edge.halfedges() {
                let hid = he.id;
                let face = he.face();
                for tid in face.data.property.tets {
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

            let first_he = edge.halfedge();
            let first_hid = first_he.id;
            let mut sorted_halfedges = Vec::with_capacity(tet_halfedges_map.len());
            sorted_halfedges.push(first_hid);
            let mut curr_hid = first_hid;
            let mut curr_tid = first_he.face().data.property.tets[0];
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

                let candidates = tet_mesh.halfedge(curr_hid).face().data.property.tets;
                curr_tid = if candidates[0] == curr_tid {
                    candidates[1]
                } else {
                    candidates[0]
                };
            }

            debug_assert!(sorted_halfedges.len() == tet_halfedges_map.len());
            for (h1, h2) in sorted_halfedges.into_iter().circular_tuple_windows() {
                unsafe { mesh_ptr.as_mut().set_he_sibling(h1, h2) };
            }
        }

        TetSet {
            mesh: tet_mesh,
            tets,
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
            if self.mesh.edge_data(eid).property.square_len != len {
                // edge changed
                continue;
            }
            split_bump.reset();

            let start_fid = self.split_edge(eid, srf_datum, &split_bump).into();
            for face in self
                .mesh
                .face_range(start_fid, self.mesh.n_faces() - *start_fid)
            {
                for tid in face.data.property.tets {
                    unsafe {
                        push_longest_edge(
                            tid,
                            NonNull::from_ref(self).as_mut(),
                            srf_datum,
                            &mut pq,
                            sq_eps,
                        );
                    }
                }
            }
        }
    }

    pub(crate) fn build_descending_vertex_links(&self) -> Vec<VertexId> {
        let mut descent_links = vec![VertexId::default(); self.mesh.n_vertices()];
        for tet_verts in self.tets.iter().map(|tet| &tet.vertices) {
            let min_idx = tet_verts
                .map(|vid| &self.mesh.vertex_data(vid).property.pt)
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
                let next_vid = descent_links[*vid];
                if !next_vid.valid() {
                    descent_links[*vid] = min_v;
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
        let n_extra_edges = self.mesh.edge(split_eid).halfedges().count();
        self.mesh.halfedge_reserve(n_extra_edges * 6);
        self.mesh.edge_reserve(n_extra_edges + 1);
        self.mesh.face_reserve(n_extra_edges * 3);

        let [va, vb] = self.mesh.edge(split_eid).vertices().map(|v| v.id);

        let new_eid = self.mesh.n_edges_capacity().into();
        let ve = self.mesh.split_edge(split_eid);
        let mut mesh: NonNull<SurfaceMesh<_, _, _, _, std::alloc::Global>> =
            NonNull::from_mut(&mut self.mesh);

        let new_pt = {
            // update points
            let mesh_mut = unsafe { mesh.as_mut() };
            let pa = &self.mesh.vertex_data(va).property.pt;
            let pb = &self.mesh.vertex_data(vb).property.pt;
            let pe = &mut mesh_mut.vertex_data_mut(ve).property.pt;
            pe[0] = (pa[0] + pb[0]) * 0.5;
            pe[1] = (pa[1] + pb[1]) * 0.5;
            pe[2] = (pa[2] + pb[2]) * 0.5;
            pe
        };
        {
            // update edge lengths
            let old_edge_data = &mut self.mesh.edge_data_mut(split_eid).property;
            old_edge_data.square_len *= 0.25;
            self.mesh.edge_data_mut(new_eid).property.square_len = old_edge_data.square_len;
        }
        let mut side_halfedges = Vec::with_capacity_in(n_extra_edges, alloc);
        side_halfedges.extend(self.mesh.edge(split_eid).halfedges().map(|he| {
            let (h_ac, h_bc, vc) = if he.to().id == vb {
                let he_next = he.next();
                let vc = he_next.data.vertex;
                (he.prev().prev().id, he_next.id, vc)
            } else {
                let he_prev = he.prev();
                let vc = he_prev.prev().data.vertex;
                (he.next().next().id, he_prev.id, vc)
            };
            let face = he.face();
            let mesh_mut = unsafe { mesh.as_mut() };
            let new_hid = mesh_mut.split_face(face.id, ve, vc);
            let mut new_he = mesh_mut.halfedge_mut(new_hid);
            let new_face = new_he.face_mut();
            new_face.data.property.tets = face.data.property.tets;
            if vc.valid() {
                new_he.edge_mut().data.property.square_len = square_norm(
                    sub_short::<3, _>(
                        self.mesh.vertex_data(vc).property.pt.as_slice(),
                        new_pt.as_slice(),
                    )
                    .as_slice(),
                );
            }
            debug_assert_eq!(new_he.data.vertex, vc);
            debug_assert_eq!(new_he.prev().data.vertex, ve);
            [new_he.id, new_he.sibling().id, h_ac, h_bc]
        }));

        let ret = self.mesh.n_faces();
        let mut prev_tid = INVALID_IND;

        let mut evaluation_map = HashMap::with_capacity_in(
            self.tets[self
                .mesh
                .edge(split_eid)
                .halfedge()
                .face()
                .data
                .property
                .tets[0]]
                .surface_evaluations
                .len(),
            alloc,
        );

        for ([h_ec, h_ce, h_ac, h_bc], [h_ed, h_de, h_ad, h_bd]) in
            side_halfedges.into_iter().circular_tuple_windows()
        {
            let left_bottom_face = unsafe { mesh.as_mut().halfedge_mut(h_ac).face_mut() };
            let mut left_top_face = unsafe { mesh.as_mut().halfedge_mut(h_bc).face_mut() };
            let right_bottom_face = unsafe { mesh.as_mut().halfedge_mut(h_ad).face_mut() };
            let mut right_top_face = unsafe { mesh.as_mut().halfedge_mut(h_bd).face_mut() };

            let [vc, vd] = [self.mesh.he_to(h_ec), self.mesh.he_to(h_ed)];

            // Optimized tetrahedron finding
            let old_tid = if prev_tid == INVALID_IND {
                let left_tets = left_bottom_face.data.property.tets;
                let right_tets = right_bottom_face.data.property.tets;
                if left_tets[0] == right_tets[0] || left_tets[0] == right_tets[1] {
                    left_tets[0]
                } else {
                    left_tets[1]
                }
            } else {
                let face_tets = left_bottom_face.data.property.tets;
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

            let bottom_hid = self.mesh.he_from_oppo_vertex(bottom_fid, va);

            let h_cd = if self.mesh.he_to(bottom_hid) != vd {
                self.mesh.halfedge(bottom_hid).twin().id
            } else {
                bottom_hid
            };

            debug_assert!(self.mesh.he_from(h_cd) == self.mesh.he_to(h_ec));
            debug_assert!(self.mesh.he_to(h_cd) == self.mesh.he_to(h_ed));

            let new_fid = unsafe {
                mesh.as_mut()
                    .add_face_by_halfedges(&[h_ec, h_cd, h_de], false)
            };
            self.mesh.face_data_mut(new_fid).property.tets = [new_tid, old_tid];

            // set sibling halfedges - optimized by caching halfedge lookups
            {
                let mut set_sibling = |prev_hid, next_hid, curr_hid| unsafe {
                    let mut prev_he = mesh.as_mut().halfedge_mut(prev_hid);
                    let mut next_he = mesh.as_mut().halfedge_mut(next_hid);
                    let mut curr_he = mesh.as_mut().halfedge_mut(curr_hid);
                    debug_assert!(curr_he.edge().id == prev_he.edge().id);
                    debug_assert!(curr_he.edge().id == next_he.edge().id);
                    if prev_he.data.sibling() == next_hid {
                        prev_he.insert_sibling(&mut curr_he);
                    } else {
                        debug_assert!(next_he.data.sibling() == prev_hid);
                        next_he.insert_sibling(&mut curr_he);
                    }
                };

                let mut new_hid = self.mesh.f_halfedge(new_fid);
                set_sibling(h_ec, h_ce, new_hid);

                new_hid = self.mesh.he_next(new_hid);
                // Optimized opposite halfedge finding
                let oppo_bottom_hid = self
                    .mesh
                    .face(top_fid)
                    .halfedges()
                    .find_map(|he| {
                        if he.next().data.vertex == vb {
                            Some(he.id)
                        } else {
                            None
                        }
                    })
                    .unwrap();
                set_sibling(bottom_hid, oppo_bottom_hid, new_hid);

                new_hid = self.mesh.he_next(new_hid);
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
                self.mesh.he_edge(h_ec), // e_ec
                self.mesh.he_edge(h_ed), // e_ed
                self.mesh.he_edge(h_cd), // e_cd
                self.mesh.he_edge(h_ac), // e_ac
                self.mesh.he_edge(h_bc), // e_bc
                self.mesh.he_edge(h_ad), // e_ad
                self.mesh.he_edge(h_bd), // e_bd
            ];

            // Batch face updates with early termination
            let update_face_tet = |face: &mut FaceMut<'_, SurfaceMesh<_, (), _, FaceData, _>>| {
                let tets = &mut face.data.property.tets;
                if tets[0] == old_tid {
                    tets[0] = new_tid;
                } else {
                    debug_assert!(tets[1] == old_tid);
                    tets[1] = new_tid;
                }
            };
            update_face_tet(&mut left_top_face);
            update_face_tet(&mut right_top_face);
            unsafe { update_face_tet(&mut mesh.as_mut().face_mut(top_fid)) };

            let new_tet = if !vc.valid() {
                let new_tet = Tet {
                    vertices: [vb, ve, vd, vc],
                    edges: [split_eid, e_bd, e_bc, e_ed, e_ec, e_cd],
                    faces: [new_fid, top_fid, left_top_face.id, right_top_face.id],
                    surface_evaluations: get_srf_eval([1, 4, 3, 2]),
                };
                tet.vertices = [ve, va, vd, vc];
                tet.edges = [new_eid, e_ed, e_ec, e_ad, e_ac, e_cd];
                tet.faces = [
                    bottom_fid,
                    new_fid,
                    left_bottom_face.id,
                    right_bottom_face.id,
                ];
                tet.surface_evaluations = get_srf_eval([4, 0, 3, 2]);
                new_tet
            } else {
                let new_tet = Tet {
                    vertices: [ve, vb, vc, vd],
                    edges: [split_eid, e_ec, e_ed, e_bc, e_bd, e_cd],
                    faces: [top_fid, new_fid, right_top_face.id, left_top_face.id],
                    surface_evaluations: get_srf_eval([4, 1, 2, 3]),
                };
                tet.vertices = [va, ve, vc, vd];
                tet.edges = [new_eid, e_ac, e_ad, e_ec, e_ed, e_cd];
                tet.faces = [
                    new_fid,
                    bottom_fid,
                    right_bottom_face.id,
                    left_bottom_face.id,
                ];
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
                        let [vc, vd] = self.mesh.e_vertices(eid);
                        debug_assert!((va == vc && vb == vd) || (va == vd && vb == vc));
                        if va.valid() && vb.valid() {
                            let pa = self.mesh.vertex_data(va).property.pt;
                            let pb = self.mesh.vertex_data(vb).property.pt;
                            let sq_len = square_norm(&sub_short::<3, _>(&pa, &pb));
                            debug_assert!(
                                (sq_len - self.mesh.edge_data(eid).property.square_len).abs()
                                    < 1e-12
                            );
                        }

                        let halfedges =
                            Vec::from_iter(self.mesh.edge(eid).halfedges().map(|he| he.id));
                        for (&h1, &h2) in halfedges.iter().circular_tuple_windows() {
                            let face1 = self.mesh.halfedge(h1).face();
                            let face2 = self.mesh.halfedge(h2).face();
                            let [t1, t2] = face1.data.property.tets;
                            let [t3, t4] = face2.data.property.tets;
                            assert!(t1 == t3 || t1 == t4 || t2 == t3 || t2 == t4);
                        }
                    }

                    if tet.vertices[3].valid() {
                        let points = tet
                            .vertices
                            .map(|vid| self.mesh.vertex_data(vid).property.pt);
                        for eval in &tet.surface_evaluations {
                            for i in 0..4 {
                                let e1 = eval.evaluation[i];
                                let e2 = srf_datum[eval.sid].surf.eval(&points[i]);
                                let length = square_norm(&sub_short::<4, _>(&e1, &e2));
                                debug_assert!(length < 1e-12);
                            }
                        }
                    }

                    let mut tet_vertices = tet.vertices.map(|vid| *vid);
                    tet_vertices.sort();
                    for (fid, vd) in tet.faces.into_iter().zip(tet.vertices) {
                        let he = self.mesh.face(fid).halfedge();
                        let mut vs0 = [
                            *he.prev().data.vertex,
                            *he.data.vertex,
                            *he.next().data.vertex,
                            *vd,
                        ];
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
    fn subdividable(&mut self, tid: usize, srf_datum: &[SurfaceData], sq_eps: f64) -> bool {
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

        let tet = &mut self.tets[tid];
        if !tet.valid()
            || tet
                .edges
                .into_iter()
                .all(|eid| self.mesh.edge_data(eid).property.square_len < sq_eps)
        {
            return false;
        }
        if tid > 1000_0000 {
            let edge_lengths = Vec::from_iter(
                tet.edges
                    .into_iter()
                    .map(|eid| self.mesh.edge_data(eid).property.square_len.sqrt()),
            );
            println!("{:?}", edge_lengths);
        }

        let tet_points = tet
            .vertices
            .map(|vid| &self.mesh.vertex_data(vid).property.pt);
        let tet_box = BBox::from_iter(tet_points);

        let mut contain_some_srf = false;
        tet.surface_evaluations.retain(|eval| {
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

        if tet.surface_evaluations.len() < 1 {
            return false;
        }

        let trans_vmat = [
            sub_short::<3, _>(tet_points[1].as_slice(), tet_points[0].as_slice()),
            sub_short::<3, _>(tet_points[2].as_slice(), tet_points[0].as_slice()),
            sub_short::<3, _>(tet_points[3].as_slice(), tet_points[0].as_slice()),
            sub_short::<3, _>(tet_points[2].as_slice(), tet_points[1].as_slice()),
            sub_short::<3, _>(tet_points[3].as_slice(), tet_points[1].as_slice()),
            sub_short::<3, _>(tet_points[3].as_slice(), tet_points[2].as_slice()),
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

        let n_surfaces = tet.surface_evaluations.len();
        let mut interpolant_vec = Vec::with_capacity(n_surfaces);
        let mut interpolant_diff_vec = Vec::with_capacity(n_surfaces);
        let mut val_diff_vec = Vec::with_capacity(n_surfaces);

        for eval in tet.surface_evaluations.iter() {
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

fn push_longest_edge(
    tid: usize,
    tets: &mut TetSet,
    srf_datum: &[SurfaceData],
    pq: &mut BinaryHeap<EdgeAndLen>,
    sq_eps: f64,
) {
    if tets.subdividable(tid, srf_datum, sq_eps) {
        let longest_eid = *tets.tets[tid]
            .edges
            .iter()
            .max_by(|&&ea, &&eb| {
                tets.mesh
                    .edge_data(ea)
                    .property
                    .square_len
                    .partial_cmp(&tets.mesh.edge_data(eb).property.square_len)
                    .unwrap()
            })
            .unwrap();

        pq.push(EdgeAndLen {
            eid: longest_eid,
            len: tets.mesh.edge_data(longest_eid).property.square_len,
        });
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
