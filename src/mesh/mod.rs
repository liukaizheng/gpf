mod element;
mod surface_mesh;

#[macro_use]
mod mesh_macro;

pub use element::*;
pub use surface_mesh::*;

use crate::INVALID_IND;

const BL_START: usize = INVALID_IND - 1;
pub enum FaceOrBoundaryLoopId {
    Face(FaceId),
    BoundaryLoop(BoundaryLoopId),
    INVALID,
}

pub trait Mesh: Sized {
    fn n_vertices(&self) -> usize;
    fn n_halfedges(&self) -> usize;
    fn n_edges(&self) -> usize;
    fn n_faces(&self) -> usize;
    fn n_boundary_loops(&self) -> usize;

    fn n_vertices_capacity(&self) -> usize;
    fn n_halfedges_capacity(&self) -> usize;
    fn n_edges_capacity(&self) -> usize;
    fn n_faces_capacity(&self) -> usize;
    fn n_boundary_loop_capacity(&self) -> usize;

    fn vertex_is_valid(&self, vid: VertexId) -> bool;
    fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool;
    fn edge_is_valid(&self, eid: EdgeId) -> bool;
    fn face_is_valid(&self, fid: FaceId) -> bool;
    fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool;

    fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self>;
    fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self>;
    fn edge<'a>(&'a self, eid: EdgeId) -> EdgeIter<'a, Self>;
    fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self>;
    fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self>;

    fn vertices(&self) -> VertexIter<'_, Self>;
    fn halfedges(&self) -> HalfedgeIter<'_, Self>;
    fn edges(&self) -> EdgeIter<'_, Self>;
    fn faces(&self) -> FaceIter<'_, Self>;

    /// the halfedge starting from this vertex
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId;

    /// the start vertex of the halfedge
    fn he_vertex(&self, hid: HalfedgeId) -> VertexId;
    /// the end vertex of the halfedge
    fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId;
    /// the next halfedge of the halfedge
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the previous halfedge of the halfedge
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the twin halfedge of the halfedge
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the sibling halfedge of the halfedge
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the next incoming halfedge of the halfedge
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the next outgoing halfedge of the halfedge
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the edge of the halfedge
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId;
    /// the face or boundary loop of the halfedge
    fn he_face_or_boundary_loop(&self, hid: HalfedgeId) -> FaceOrBoundaryLoopId;

    /// the first halfedge of the edge
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId;

    /// the first halfedge of the face
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId;

    fn use_implicit_twin(&self) -> bool;

    fn new_halfedge(&mut self, is_interior: bool) -> HalfedgeId;

    fn new_edge(&mut self) -> EdgeId;
}
