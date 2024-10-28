mod element;
mod manifold_mesh;
mod mesh;
mod mesh_core_data;
mod surface_mesh;

use std::alloc::Allocator;

pub use element::*;
pub use manifold_mesh::*;
pub use mesh::*;
pub use surface_mesh::*;

use crate::point;

#[inline]
pub fn square_edge_length<M: Mesh>(points: &[f64], eid: EdgeId, mesh: &M) -> f64 {
    let [va, vb] = mesh.e_vertices(eid);
    let pa = point(points, va.0);
    let pb = point(points, vb.0);
    pa.iter().zip(pb).map(|(a, b)| a - b).map(|x| x * x).sum()
}

#[inline]
pub(crate) fn clone_vec_in<T: Clone, A: Allocator + Copy>(vec: &[T], alloc: A) -> Vec<T, A> {
    let mut new_vec = Vec::with_capacity_in(vec.len(), alloc);
    new_vec.extend_from_slice(vec);
    new_vec
}