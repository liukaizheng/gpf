use std::alloc::Allocator;

use crate::{
    mesh::{square_edge_length, EdgeId, FaceId, Mesh, SurfaceMesh, VertexId},
    point, INVALID_IND,
};

pub(crate) struct TetSet {
    pub(crate) mesh: SurfaceMesh,
    pub(crate) tet_vertices: Vec<[VertexId; 4]>,
    pub(crate) tet_edges: Vec<[EdgeId; 6]>,
    pub(crate) tet_faces: Vec<[FaceId; 4]>,
    pub(crate) face_tets: Vec<[usize; 2]>,
    pub(crate) points: Vec<f64>,
    pub(crate) square_edge_lengths: Vec<f64>,
}

impl TetSet {
    #[inline]
    pub(crate) fn vertices_in<A: Allocator + Copy>(
        &self,
        tid: usize,
        alloc: A,
    ) -> Vec<VertexId, A> {
        self.tet_vertices[tid].to_vec_in(alloc)
    }

    #[inline]
    pub(crate) fn tet_face_reversed(&self, tid: usize, fid: FaceId) -> bool {
        let tets = self.face_tets[fid];
        debug_assert!(tets[0] == tid || tets[1] == tid);
        tid == tets[1]
    }

    pub(crate) fn split_edge<A: Allocator + Copy>(
        &mut self,
        eid: EdgeId,
        alloc: A,
    ) -> (VertexId, Vec<usize, A>) {
        let [va, vb] = self.mesh.e_vertices(eid);
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
        let mut oppo_halfedges = Vec::new_in(alloc);
        let mut bottom_faces = Vec::new_in(alloc);
        let mut top_faces = Vec::new_in(alloc);
        for (&tid, &[fa, fb]) in &tet_faces_map {
            let bottom_top_faces = || {
                let mut fc = FaceId::default();
                let mut fd = FaceId::default();
                for fid in self.tet_faces[tid] {
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
            if self.tet_face_reversed(tid, bottom_fid) {
                oppo_halfedges.push(self.mesh.he_twin(self.mesh.he_next(hid)));
            } else {
                oppo_halfedges.push(self.mesh.he_next(hid));
            }
        }
        let new_vert = self.mesh.split_edge(eid, alloc);
        let [bottom_eid, top_eid] = {
            let hid = self.mesh.e_halfedge(eid);
            if self.mesh.he_vertex(hid) == va || self.mesh.he_tip_vertex(hid) == va {
                [eid, self.square_edge_lengths.len().into()]
            } else {
                [self.square_edge_lengths.len().into(), eid]
            }
        };
        {
            // update points
            let pa = point(&self.points, va.0);
            let pb = point(&self.points, vb.0);
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

        let mut new_halfedges = Vec::new_in(alloc);
        new_halfedges.extend(oppo_verts.into_iter().zip(&faces).map(|(v, &fid)| {
            let hid = self.mesh.split_face(fid, new_vert, v, alloc);
            self.square_edge_lengths.push(square_edge_length(
                &self.points,
                self.mesh.he_edge(hid),
                &self.mesh,
            ));
            self.face_tets.push(self.face_tets[fid]);
            if self.mesh.he_vertex(hid) == new_vert {
                hid
            } else {
                self.mesh.he_twin(hid)
            }
        }));

        debug_assert!(new_halfedges
            .iter()
            .all(|&hid| self.mesh.he_vertex(hid) == new_vert));

        let mut result_tets = Vec::with_capacity_in(tet_faces_map.len() << 1, alloc);
        for ((((tid, face_indices), oppo_hid), bottom_fid), top_fid) in tet_faces_map
            .into_iter()
            .zip(oppo_halfedges)
            .zip(bottom_faces)
            .zip(top_faces)
        {
            let mut ha = new_halfedges[face_indices[0]];
            let mut hb = new_halfedges[face_indices[1]];
            if self.mesh.he_tip_vertex(ha) != self.mesh.he_vertex(oppo_hid) {
                std::mem::swap(&mut ha, &mut hb);
            }
            hb = self.mesh.he_twin(hb);

            debug_assert!(self.mesh.he_tip_vertex(ha) == self.mesh.he_vertex(oppo_hid));
            debug_assert!(self.mesh.he_vertex(hb) == self.mesh.he_tip_vertex(oppo_hid));
            debug_assert!(self.mesh.he_tip_vertex(hb) == self.mesh.he_vertex(ha));

            let new_fid = self.mesh.add_face_by_halfedges(&[ha, oppo_hid, hb], alloc);
            let new_tid = self.tet_vertices.len();
            self.face_tets.push([tid, new_tid]);

            let mut ot_faces = Vec::with_capacity_in(4, alloc); // old tet faces
            let mut nt_faces = Vec::with_capacity_in(4, alloc); // new tet faces

            let mut ot_edges = Vec::with_capacity_in(6, alloc); // old tet edges
            let mut nt_edges = Vec::with_capacity_in(6, alloc); // new tet edges
            ot_edges.push(top_eid);
            nt_edges.push(bottom_eid);

            let replace = |tets: &mut [usize; 2]| {
                for t in tets {
                    if *t == tid {
                        *t = new_tid;
                        break;
                    }
                }
            };

            // bottom face
            nt_faces.push(bottom_fid);
            replace(&mut self.face_tets[bottom_fid]);

            // old tet top face
            ot_faces.push(top_fid);

            // new face
            ot_faces.push(new_fid);
            nt_faces.push(new_fid);

            // fa and fb
            for idx in face_indices {
                let fid = faces[idx];
                // when the face is split, the f_halfedge of the face is the split halfedge
                let he = self.mesh.face(fid).halfedge();
                if *he.next().to() == va {
                    debug_assert!(*he.twin().next().to() == vb);
                    nt_faces.push(fid);
                    replace(&mut self.face_tets[fid]);
                    ot_faces.push(self.mesh.he_face(*he.twin()));
                } else {
                    debug_assert!(*he.next().to() == vb);
                    debug_assert!(*he.twin().next().to() == va);
                    let twin_fid = self.mesh.he_face(*he.twin());
                    nt_faces.push(twin_fid);
                    replace(&mut self.face_tets[twin_fid]);
                    ot_faces.push(fid);
                }
            }
            debug_assert!(ot_faces.len() == 4);
            debug_assert!(nt_faces.len() == 4);

            let collect_edges = |vec: &mut Vec<EdgeId, _>, candidate_faces, exclude_eid| {
                for &fid in candidate_faces {
                    for he in self.mesh.face(fid).halfedges() {
                        let e = he.edge();
                        if *e != exclude_eid {
                            vec.push(*e);
                        }
                    }
                }
            };
            collect_edges(&mut ot_edges, &ot_faces[2..], top_eid);
            collect_edges(&mut nt_edges, &nt_faces[2..], bottom_eid);

            let oppo_eid = self.mesh.he_edge(oppo_hid);
            ot_edges.push(oppo_eid);
            nt_edges.push(oppo_eid);

            debug_assert!(ot_edges.len() == 6);
            debug_assert!(nt_edges.len() == 6);

            let [vc, vd] = self.mesh.e_vertices(self.mesh.he_edge(oppo_hid));

            self.tet_vertices[tid] = [new_vert, vc, vd, vb];
            self.tet_edges[tid] = [
                ot_edges[0],
                ot_edges[1],
                ot_edges[2],
                ot_edges[3],
                ot_edges[4],
                ot_edges[5],
            ];
            self.tet_faces[tid] = [ot_faces[0], ot_faces[1], ot_faces[2], ot_faces[3]];

            debug_assert!(ot_faces.iter().all(|&fid| {
                let tets = self.face_tets[fid];
                tets[0] == tid || tets[1] == tid
            }));
            debug_assert!(nt_faces.iter().all(|&fid| {
                let tets = self.face_tets[fid];
                tets[0] == new_tid || tets[1] == new_tid
            }));

            self.tet_vertices.push([new_vert, va, vc, vd]);
            self.tet_edges.push([
                nt_edges[0],
                nt_edges[1],
                nt_edges[2],
                nt_edges[3],
                nt_edges[4],
                nt_edges[5],
            ]);
            self.tet_faces
                .push([nt_faces[0], nt_faces[1], nt_faces[2], nt_faces[3]]);
            result_tets.extend([tid, new_tid]);
        }
        (new_vert, result_tets)
    }
}
