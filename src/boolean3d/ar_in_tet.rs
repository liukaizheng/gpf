use std::{
    alloc::{Allocator, Global},
    cell::LazyCell,
};

use crate::{
    mesh::{clone_vec_in, ElementId, FaceId, Mesh, SurfaceMesh, VertexId},
    predicates::{det4, sign_reverse, sign_reversed, Orientation}, INVALID_IND,
};

#[derive(Clone)]
struct FaceData {
    /// plane id
     pid: usize,
    /// `cells[0]`: inner cell of this face
    /// `cells[1]`: outer cell of this face
    cells: [usize; 2],
}

pub(super) struct Arrangement<A: Allocator + Copy> {
    mesh: SurfaceMesh<A>,
    vertices: Vec<[usize; 3], A>,
    edges: Vec<[usize; 2], A>,
    planes: Vec<[f64; 4], A>,
    face_data: Vec<FaceData, A>,
    cell_faces: Vec<Vec<FaceId, A>, A>,
}

fn orient3d<A: Allocator + Copy>(
    pa: &[f64; 4],
    pb: &[f64; 4],
    pc: &[f64; 4],
    pd: &[f64; 4],
    alloc: A,
) -> Orientation {

    #[rustfmt::skip]
    let numerator_sign = det4(
        pa[0], pa[1], pa[2], pa[3],
        pb[0], pb[1], pb[2], pb[3],
        pc[0], pc[1], pc[2], pc[3],
        pd[0], pd[1], pd[2], pd[3],
        alloc,
    );
    let denominator_sign = det4(
        pa[0], pa[1], pa[2], pa[3],
        pb[0], pb[1], pb[2], pb[3],
        pc[0], pc[1], pc[2], pc[3],
        1.0  , 1.0  , 1.0  ,   1.0,
        alloc,
    );

    if denominator_sign == Orientation::Positive {
        return numerator_sign;
    } else {
        debug_assert!(denominator_sign != Orientation::Zero);
        return sign_reverse(numerator_sign);
    }
}

impl<A: Allocator + Copy> Arrangement<A> {
    fn new_tet(alloc: A) -> Self {
        let mesh = SurfaceMesh::new([[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]], alloc);
        let mut vertices = Vec::with_capacity_in(4, alloc);
        vertices.extend(mesh.vertices().map(|v| {
            let mut res = [0; 3];
            for (fid, hid) in res.iter_mut().zip(v.incoming_halfedges()) {
                *fid = mesh.he_face(*hid).0;
            }
            res
        }));

        let mut edges = Vec::with_capacity_in(6, alloc);
        edges.extend(mesh.edges().map(|e| {
            let mut res = [0; 2];
            for (fid, hid) in res.iter_mut().zip(e.halfedges()) {
                *fid = mesh.he_face(*hid).0;
            }
            res
        }));

        let mut planes = Vec::with_capacity_in(4, alloc);
        planes.push([1.0, 0.0, 0.0, 0.0]);
        planes.push([0.0, 1.0, 0.0, 0.0]);
        planes.push([0.0, 0.0, 1.0, 0.0]);
        planes.push([0.0, 0.0, 0.0, 1.0]);

        let mut face_data = Vec::with_capacity_in(4, alloc);
        face_data.extend((0..4).map(|pid| FaceData {
            pid,
            cells: [0, INVALID_IND],
        }));

        let mut one_cell_faces = Vec::with_capacity_in(4, alloc);
        one_cell_faces.extend((0..4).map(|idx| FaceId::from(idx)));

        let mut cell_faces = Vec::with_capacity_in(1, alloc);
        cell_faces.push(one_cell_faces);

        Self {
            mesh,
            vertices,
            edges,
            planes,
            face_data,
            cell_faces,
        }
    }

    fn clone_in<A1: Allocator + Copy>(&self, alloc: A1) -> Arrangement<A1> {
        let mut cell_faces = Vec::with_capacity_in(self.cell_faces.len(), alloc);
        for v in &self.cell_faces {
            cell_faces.push(clone_vec_in(v, alloc));
        }
        Arrangement {
            mesh: self.mesh.clone_in(alloc),
            vertices: clone_vec_in(&self.vertices, alloc),
            edges: clone_vec_in(&self.edges, alloc),
            planes: clone_vec_in(&self.planes, alloc),
            face_data: clone_vec_in(&self.face_data, alloc),
            cell_faces,
        }
    }

    #[inline]
    fn is_face_inner_cell(&self, fid: FaceId, cid: usize) -> bool {
        let cells = &self.face_data[fid].cells;
        debug_assert!(cells[0] == cid || cells[1] == cid);
        cells[0] == cid
    }

    fn add_plane<A1: Allocator + Copy>(&mut self, plane: &[f64; 4], alloc: A1) {
        let pid = self.planes.len();
        self.planes.push(plane.clone());
        let mut vert_orientations = Vec::with_capacity_in(self.mesh.n_vertices_capacity(), alloc);
        vert_orientations.extend(self.mesh.vertices().map(|v| {
            let vid = *v;
            let v_planes = &self.vertices[vid];
            orient3d(&self.planes[v_planes[0]], &self.planes[v_planes[1]],  &self.planes[v_planes[2]], plane, alloc)
        }));

        for f in self.mesh.faces() {
            if f.vertices().all(|v| vert_orientations[*v] == Orientation::Zero) {
                println!("add same plane");
                return;
            }
        }

        self.split_edges(&mut vert_orientations, pid, alloc);
        self.split_faces(&vert_orientations, pid);
        self.split_cells(&vert_orientations, pid, alloc);
    }


    fn split_edges<A1: Allocator + Copy>(&mut self, vert_orientations: &mut Vec<Orientation, A1>, pid: usize, alloc: A1) {
        let n_old_edges = self.mesh.n_edges_capacity();
        for eid in 0..n_old_edges {
            let eid = eid.into();
            let [va, vb] = self.mesh.e_vertices(eid);
            if sign_reversed(vert_orientations[va], vert_orientations[vb]) {
                self.mesh.split_edge(eid, alloc);
                let e_planes = &self.edges[eid];
                self.vertices.push([e_planes[0], e_planes[1], pid]);
                vert_orientations.push(Orientation::Zero);
                self.edges.push(e_planes.clone());
            }
        }
    }

    fn split_faces(&mut self, vert_orientations: &[Orientation], pid: usize) {
        let n_old_faces = self.mesh.n_faces_capacity();
        for fid in 0..n_old_faces {
            let fid = fid.into();
            let mut first_zero_vid = VertexId::default();
            let first_hid = self.mesh.f_halfedge(fid);
            let mut curr_hid = first_hid;
            let mut prev_ori = vert_orientations[self.mesh.he_from(curr_hid)];
            if prev_ori == Orientation::Zero {
                first_zero_vid = self.mesh.he_from(curr_hid);
            }
            loop {
                let vid = self.mesh.he_to(curr_hid);
                let ori = vert_orientations[vid];
                let next_hid = self.mesh.he_next(curr_hid);
                if ori == Orientation::Zero {
                    if prev_ori == Orientation::Zero {
                        self.mesh.set_f_halfedge(fid, curr_hid);
                        break;
                    } else {
                        if first_zero_vid.valid() {
                            if self.mesh.he_to(next_hid) == first_zero_vid {
                                self.mesh.set_f_halfedge(fid, next_hid);
                                break;
                            } else if vid == first_zero_vid {
                                break;
                            }
                            self.mesh.split_face(fid, first_zero_vid, vid);
                            let new_fid = self.face_data.len().into();
                            let data = self.face_data[fid].clone();
                            self.edges.push([pid, data.pid]);
                            for cid in data.cells {
                                if cid != INVALID_IND {
                                    self.cell_faces[cid].push(new_fid);
                                }
                            }
                            self.face_data.push(data);
                            break;
                        } else {
                            first_zero_vid = vid;
                        }
                    }
                }
                curr_hid = next_hid;
                if curr_hid == first_hid {
                    break;
                }

                prev_ori = ori;
            }
        }
    }

    fn split_cells<A1: Allocator + Copy>(&mut self, vert_orientations: &[Orientation], pid: usize, alloc: A1) {
        let n_old_cells = self.cell_faces.len();
        let self_alloc = self.vertices.allocator().clone();
        for cid in 0..n_old_cells {
            let mut pos_cell_faces = Vec::new_in(self_alloc);
            let mut neg_cell_faces = Vec::new_in(self_alloc);

            let mut start_zero_vid = VertexId::default();
            let mut n_halfedges = 0;
            for &fid in &self.cell_faces[cid] {
                let hid = self.mesh.f_halfedge(fid);
                let [va, vb] = self.mesh.he_vertices(hid);
                let mut non_zero_vid = VertexId::default();
                if vert_orientations[va] == Orientation::Zero {
                    if vert_orientations[vb] == Orientation::Zero {
                        let vc = *self.mesh.halfedge(hid).next().to();
                        debug_assert!(vert_orientations[vc] != Orientation::Zero);

                        let is_pos = vert_orientations[vc] == Orientation::Positive;
                        if is_pos {
                            pos_cell_faces.push(fid);
                        } else {
                            neg_cell_faces.push(fid);
                        }

                        if is_pos == self.is_face_inner_cell(fid, cid) {
                            self.mesh.set_v_halfedge(va, hid);
                            if !start_zero_vid.valid() {
                                start_zero_vid = va;
                            }
                            n_halfedges += 1;
                        }
                    } else {
                        non_zero_vid = vb;
                    }
                } else {
                    non_zero_vid = va;
                }

                if non_zero_vid.valid() {
                    if vert_orientations[non_zero_vid] == Orientation::Positive {
                        pos_cell_faces.push(fid);
                    } else {
                        neg_cell_faces.push(fid);
                    }
                }
            }

            if !start_zero_vid.valid() {
                continue;
            }

            let mut curr_vid = start_zero_vid;
            let mut new_halfedges = Vec::with_capacity_in(n_halfedges, alloc);
            loop {
                let hid = self.mesh.v_halfedge(curr_vid);
                new_halfedges.push(hid);
                curr_vid = self.mesh.he_to(hid);
                if curr_vid == start_zero_vid {
                    break;
                }
            }

            let new_fid = self.mesh.add_face_by_halfedges(&new_halfedges);
            let new_cid = self.cell_faces.len();
            self.face_data.push(FaceData { pid , cells: [cid, new_cid] });

            for &fid in &pos_cell_faces {
                for c in &mut self.face_data[fid].cells {
                    if *c == cid {
                        *c = new_cid;
                    }
                }
            }

            neg_cell_faces.push(new_fid);
            pos_cell_faces.push(new_fid);
            self.cell_faces[cid] = neg_cell_faces;
            self.cell_faces.push(pos_cell_faces);
        }
    }

}

pub(super) fn arrangement_for_tet<A: Allocator + Copy>(planes: &[[f64; 4]], alloc: A) {
    let one_tet = LazyCell::new(|| Arrangement::new_tet(Global));
    let mut ar = one_tet.clone_in(Global);
    for plane in planes {
        ar.add_plane(plane, alloc);
    }
}
