use std::alloc::Allocator;

use hashbrown::HashMap;
use tinyvec::TinyVec;

use crate::{
    INVALID_IND,
    geometry::{BBox, Surf},
    mesh::{
        EdgeId, ElementId, FaceId, HalfedgeId, Mesh, SurfaceMesh, VertexId, square_edge_length,
    },
    point_3,
};

#[derive(Default)]
pub(crate) struct SurfaceEvaluation {
    pub(crate) sid: usize,
    pub(crate) evaluation: [[f64; 4]; 4],
}

pub(crate) struct Tet {
    pub(crate) vertices: [VertexId; 4],
    pub(crate) edges: [HalfedgeId; 6],
    pub(crate) faces: [FaceId; 4],
    pub(crate) surface_evaluations: TinyVec<[SurfaceEvaluation; 3]>,
}

pub(crate) struct TetSet {
    pub(crate) mesh: SurfaceMesh<std::alloc::Global>,
    pub(crate) tets: Vec<Tet>,
    pub(crate) face_tets: Vec<[usize; 2]>,
    pub(crate) points: Vec<f64>,
    pub(crate) square_edge_lengths: Vec<f64>,
    pub(crate) surf_indices: Vec<TinyVec<[usize; 3]>>,
}

#[inline]
pub(crate) fn tet_face_reversed(face_tets: &[usize; 2], tid: usize) -> bool {
    debug_assert!(face_tets[0] == tid || face_tets[1] == tid);
    tid == face_tets[1]
}

impl TetSet {
    pub(crate) fn form_bbox(bbox: BBox, n_surfaces: usize) {
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

        let tet_mesh = SurfaceMesh::new(
            TETS.map(|tet_vertices| TET_FACES.map(|face| face.map(|idx| tet_vertices[idx])))
                .iter()
                .flatten(),
            std::alloc::Global,
        );

        println!("the tet face is {}", tet_mesh.n_faces());

        // let hash_tri = |mut verts: [usize; 3]| {
        //     verts.sort();
        //     if verts[2] == INVALID_IND {
        //         verts[2] = 8;
        //     }
        //     (verts[2] << 6) | (verts[1] << 3) | verts[0]
        // };

        /*let mut edge_square_len_map = HashMap::<usize, f64>::with_capacity(19);
        let mut face_tet_map = HashMap::<usize, [TetHandle; 2]>::with_capacity(18);

        let mut tets = Vec::from_iter(TETS.into_iter().enumerate().map(|(tid, tet_verts)| {
            let mut edge_square_lengths = [0.0; 6];
            for ((va, vb), edge_len) in tet_verts
                .into_iter()
                .tuple_combinations()
                .zip(&mut edge_square_lengths)
            {
                *edge_len = *edge_square_len_map
                    .entry(hash_seg(va, vb))
                    .or_insert(square_norm(&sub_short::<3, _>(
                        point::<3>(&points, va),
                        point::<3>(&points, vb),
                    )));
            }
            let tet_points = tet_verts.map(|idx| point::<3>(&points, idx));
            let surface_evaluations = TinyVec::from_iter(surfaces.iter().enumerate().map(
                |(sid, srf)| SurfaceEvaluation {
                    sid,
                    evaluation: tet_points.map(|p| srf.eval(p)),
                },
            ));
            for (i, face_verts) in TET_FACES
                .map(|face_verts| face_verts.map(|idx| tet_verts[idx]))
                .into_iter()
                .enumerate()
            {
                match face_tet_map.entry(hash_tri(face_verts)) {
                    hashbrown::hash_map::Entry::Occupied(mut entry) => {
                        let handle = &mut entry.get_mut()[1];
                        handle.ver = i;
                        handle.tid = tid;
                    }
                    hashbrown::hash_map::Entry::Vacant(entry) => {
                        entry.insert([TetHandle { tid, ver: i }, TetHandle::default()]);
                    }
                }
            }

            Tet {
                vertices: tet_verts,
                neighbors: Default::default(),
                edge_square_lengths,
                surface_evaluations,
            }
        }));

        let mut bdy_face_tet_map = HashMap::<usize, [TetHandle; 2]>::new();

        for [h1, h2] in face_tet_map.values_mut() {
            if h2.tid != INVALID_IND {
                continue;
            }
            let new_tid = tets.len();
            h2.tid = new_tid;
            h2.ver = 3;
            let vertices = {
                let tet = &tets[h1.tid];
                let vs = TET_FACES[h1.ver].map(|v| tet.vertices[v]);
                [vs[0], vs[2], vs[1], INVALID_IND]
            };

            for (ver, face_verts) in TET_FACES[..3].iter().enumerate() {
                let face_verts = face_verts.map(|v| vertices[v]);
                match bdy_face_tet_map.entry(hash_tri(face_verts)) {
                    hashbrown::hash_map::Entry::Occupied(mut entry) => {
                        let handle = &mut entry.get_mut()[1];
                        handle.ver = ver;
                        handle.tid = new_tid;
                    }
                    hashbrown::hash_map::Entry::Vacant(entry) => {
                        entry.insert([TetHandle { tid: new_tid, ver }, TetHandle::default()]);
                    }
                }
            }
            tets.push(Tet {
                vertices,
                neighbors: Default::default(),
                edge_square_lengths: [0.0; 6],
                surface_evaluations: TinyVec::new(),
            });
        }

        face_tet_map.extend(bdy_face_tet_map);

        for [h1, h2] in face_tet_map.into_values() {
            bond_face(&mut tets, h1, h2);
        }

        TetComplex { points, tets }*/
    }
    pub(crate) fn tet_edge_index(pa: usize, pb: usize) -> usize {
        let min_idx = if pa < pb { pa } else { pb };
        (if min_idx != 0 { 0 } else { 1 }) + 5 - pa - pb
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
        eid: EdgeId,
        alloc: A,
    ) -> (VertexId, Vec<[usize; 2], A>) {
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
