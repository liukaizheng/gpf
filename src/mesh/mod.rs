mod element;
mod manifold_mesh;
mod mesh;
mod mesh_core_data;
mod surface_mesh;

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
