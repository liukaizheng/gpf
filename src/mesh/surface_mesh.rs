use crate::{mesh::Mesh, build_connect_info, INVALID_IND};

use super::{VertexId, HalfedgeId, FaceId, BoundaryLoopId, VertexIter, HalfedgeIter, FaceIter, BoundaryLoopIter};

struct SurfaceMesh {
	he_next_arr: Vec<usize>,
	he_vertex_arr: Vec<usize>,
	he_face_arr: Vec<usize>,
	v_halfedge_arr: Vec<usize>,
	f_halfedge_arr: Vec<usize>,
	bl_halfedge_arr: Vec<usize>,

	n_vertices: usize,
	n_halfedges: usize,
	n_faces: usize,
	n_loops: usize,
}

impl Mesh for SurfaceMesh {
	build_connect_info!();

	fn n_halfedges(&self) -> usize {
		return self.n_halfedges;
	}

	fn n_faces(&self) -> usize {
		return self.n_faces;
	}

	fn n_boundary_loop(&self) -> usize {
		return self.n_loops;
	}

	fn n_vertices_capacity(&self) -> usize {
		return self.v_halfedge_arr.len();
	}

	fn n_halfedges_capacity(&self) -> usize {
		return self.he_next_arr.len();
	}

	fn n_faces_capacity(&self) -> usize {
		return self.f_halfedge_arr.len();
	}

	fn n_boundary_loop_capacity(&self) -> usize {
		return self.bl_halfedge_arr.len();
	}

	fn vertex_is_valid(&self, vid: VertexId) -> bool {
		return self.v_halfedge_arr[vid] != INVALID_IND;
	}

	fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool {
		return self.he_next_arr[hid] != INVALID_IND;
	}

	fn face_is_valid(&self, fid: FaceId) -> bool {
		return self.f_halfedge_arr[fid] != INVALID_IND;
	}

	fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool {
		return self.bl_halfedge_arr[blid] != INVALID_IND;
	}

	fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self> {
		VertexIter::new(vid, self)
	}

	fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self> {
		HalfedgeIter::new(hid, self)
	}

	fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self> {
		FaceIter::new(fid, self)
	}

	fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self> {
		BoundaryLoopIter::new(blid, self)
	}

}