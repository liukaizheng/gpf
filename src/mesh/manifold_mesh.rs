use std::{cell::RefCell, rc::Weak};

use crate::build_connect_info;

use super::{
    BoundaryLoopId, BoundaryLoopIter, EdgeData, EdgeId, EdgeIter, ElementId, FaceData, FaceId,
    FaceIter, FaceOrBoundaryLoopId, HalfedgeData, HalfedgeId, HalfedgeIter, Mesh, MeshData,
    VertexData, VertexId, VertexIter,
};

use bumpalo::Bump;

pub struct ManifoldMesh {
    v_halfedge_arr: Vec<HalfedgeId>,
    he_next_arr: Vec<HalfedgeId>,
    he_vertex_arr: Vec<VertexId>,
    he_face_arr: Vec<FaceOrBoundaryLoopId>,
    f_halfedge_arr: Vec<HalfedgeId>,
    bl_halfedge_arr: Vec<BoundaryLoopId>,

    n_vertices: usize,
    n_halfedges: usize,
    n_interior_halfedges: usize,
    n_faces: usize,
    n_boundary_loops: usize,

    pub vertices_data: Vec<Weak<RefCell<dyn MeshData<Id = VertexId>>>>,
    pub halfedges_data: Vec<Weak<RefCell<dyn MeshData<Id = HalfedgeId>>>>,
    pub edges_data: Vec<Weak<RefCell<dyn MeshData<Id = EdgeId>>>>,
    pub faces_data: Vec<Weak<RefCell<dyn MeshData<Id = FaceId>>>>,
}

impl ManifoldMesh {
    #[inline]
    pub fn n_interior_halfedges(&self) -> usize {
        self.n_interior_halfedges
    }
    /// the number of boundary loops
    #[inline(always)]
    pub fn n_boundary_loops(&self) -> usize {
        return self.n_boundary_loops;
    }

    /// the capacity of boundary loops
    #[inline(always)]
    pub fn n_boundary_loop_capacity(&self) -> usize {
        return self.bl_halfedge_arr.len();
    }

    /// boundary loop id is valid
    #[inline(always)]
    pub fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool {
        self.bl_halfedge_arr[blid].valid()
    }

    /// start boundary loop iterator
    #[inline(always)]
    pub fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a> {
        BoundaryLoopIter::new(blid, self)
    }
}

impl Mesh for ManifoldMesh {
    build_connect_info!();

    /// the number of edges
    #[inline(always)]
    fn n_edges(&self) -> usize {
        self.n_halfedges >> 1
    }

    /// the capacity of edges
    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        self.n_halfedges_capacity() >> 1
    }

    /// the edge is valid
    #[inline(always)]
    fn edge_is_valid(&self, eid: EdgeId) -> bool {
        self.he_next_arr[self.e_halfedge(eid)].valid()
    }

    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        true
    }

    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        (*hid >> 1).into()
    }

    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        (*hid ^ 1).into()
    }

    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_twin(hid)
    }

    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_twin(self.he_next(hid))
    }

    #[inline(always)]
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_next(self.he_twin(hid))
    }

    #[inline(always)]
    fn he_face_or_boundary_loop(&self, hid: HalfedgeId) -> FaceOrBoundaryLoopId {
        self.he_face_arr[hid]
    }

    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        (*eid << 1).into()
    }
}
