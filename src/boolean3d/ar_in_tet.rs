use std::{
    alloc::{Allocator, Global},
    cell::LazyCell,
};

use crate::mesh::{clone_vec_in, FaceId, Mesh, SurfaceMesh};

pub(super) struct Arrangement<A: Allocator + Copy> {
    mesh: SurfaceMesh<A>,
    vertices: Vec<[FaceId; 3], A>,
    edges: Vec<[FaceId; 2], A>,
    planes: Vec<[f64; 4], A>,
}

impl<A: Allocator + Copy> Arrangement<A> {
    fn new_tet(alloc: A) -> Self {
        let mesh = SurfaceMesh::new([[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]], alloc);
        let mut vertices = Vec::with_capacity_in(4, alloc);
        vertices.extend(mesh.vertices().map(|v| {
            let mut res = [FaceId::default(); 3];
            for (fid, hid) in res.iter_mut().zip(v.incoming_halfedges()) {
                *fid = mesh.he_face(*hid);
            }
            res
        }));

        let mut edges = Vec::with_capacity_in(6, alloc);
        edges.extend(mesh.edges().map(|e| {
            let mut res = [FaceId::default(); 2];
            for (fid, hid) in res.iter_mut().zip(e.halfedges()) {
                *fid = mesh.he_face(*hid);
            }
            res
        }));

        let mut planes = Vec::with_capacity_in(3, alloc);
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
}

pub(super) fn arrangement_for_tet<A: Allocator + Copy>(planes: &[[f64; 4]], alloc: A) {
    let one_tet = LazyCell::new(|| Arrangement::new_tet(Global));
    let arrangement = one_tet.clone_in(alloc);
}
