use std::{
    alloc::{Allocator, Global},
    cell::LazyCell,
};

use crate::{
    mesh::{clone_vec_in, FaceId, Mesh, SurfaceMesh},
    predicates::{det4, sign_reverse, sign_reversed, Orientation},
};

pub(super) struct Arrangement<A: Allocator + Copy> {
    mesh: SurfaceMesh<A>,
    vertices: Vec<[usize; 3], A>,
    edges: Vec<[usize; 2], A>,
    planes: Vec<[f64; 4], A>,
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

        Self {
            mesh,
            vertices,
            edges,
            planes,
        }
    }

    fn clone_in<A1: Allocator + Copy>(&self, alloc: A1) -> Arrangement<A1> {
        Arrangement {
            mesh: self.mesh.clone_in(alloc),
            vertices: clone_vec_in(&self.vertices, alloc),
            edges: clone_vec_in(&self.edges, alloc),
            planes: clone_vec_in(&self.planes, alloc),
        }
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

        let n_old_edges = self.mesh.n_edges_capacity();
        for eid in 0..n_old_edges {
            let eid = eid.into();
            let [va, vb] = self.mesh.e_vertices(eid);
            if sign_reversed(vert_orientations[va], vert_orientations[vb]) {
                self.mesh.split_edge(eid, alloc);
                let e_planes = &self.edges[eid];
                self.vertices.push([e_planes[0], e_planes[1], pid]);
                self.edges.push(e_planes.clone());
            }
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
