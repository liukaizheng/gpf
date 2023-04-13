mod element;
mod surface_mesh;

#[macro_use]
mod mesh_macro;

pub use element::*;
pub trait Mesh: Sized {
	fn n_vertices(&self) -> usize;
	fn n_halfedges(&self) -> usize;
	fn n_faces(&self) -> usize;
	fn n_boundary_loop(&self) -> usize;

	fn n_vertices_capacity(&self) -> usize;
	fn n_halfedges_capacity(&self) -> usize;
	fn n_faces_capacity(&self) -> usize;
	fn n_boundary_loop_capacity(&self) -> usize;

	fn vertex_is_valid(&self, vid: Vertex) -> bool;
	fn halfedge_is_valid(&self, hid: Halfedge) -> bool;
	fn face_is_valid(&self, fid: Face) -> bool;
	fn boundary_loop_is_valid(&self, blid: BoundaryLoop) -> bool;
}
