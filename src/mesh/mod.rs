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

	fn vertex_is_valid(&self, vid: VertexId) -> bool;
	fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool;
	fn face_is_valid(&self, fid: FaceId) -> bool;
	fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool;

	fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self>;
	fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self>;
	fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self>;
	fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self>;
}
