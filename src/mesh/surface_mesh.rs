use crate::{build_connect_info, mesh::Mesh, INVALID_IND};

use super::{
    BoundaryLoopId, BoundaryLoopIter, EdgeId, EdgeIter, FaceId, FaceIter, FaceOrBoundaryLoopId,
    HalfedgeId, HalfedgeIter, VertexId, VertexIter, BL_START,
};

struct SurfaceMesh {
    v_halfedge_arr: Vec<usize>,
    he_next_arr: Vec<usize>,
    he_vertex_arr: Vec<usize>,
    he_face_arr: Vec<usize>,
    f_halfedge_arr: Vec<usize>,
    bl_halfedge_arr: Vec<usize>,

    n_vertices: usize,
    n_halfedges: usize,
    n_edges: usize,
    n_faces: usize,
    n_loops: usize,

    v_he_in_start_arr: Vec<usize>,
    v_he_out_start_arr: Vec<usize>,
    he_edge_arr: Vec<usize>,
    he_vert_in_next_arr: Vec<usize>,
    he_vert_in_prev_arr: Vec<usize>,
    he_vert_out_next_arr: Vec<usize>,
    he_vert_out_prev_arr: Vec<usize>,
    he_sibling_arr: Vec<usize>,
    e_halfedge_arr: Vec<usize>,
}

impl Mesh for SurfaceMesh {
    build_connect_info!();

    /// the number of halfedges
    #[inline(always)]
    fn n_halfedges(&self) -> usize {
        return self.n_halfedges;
    }

    /// the number of edges
    #[inline(always)]
    fn n_edges(&self) -> usize {
        return self.n_edges;
    }

    /// the number of faces
    #[inline(always)]
    fn n_faces(&self) -> usize {
        return self.n_faces;
    }

    /// the number of boundary loops
    #[inline(always)]
    fn n_boundary_loops(&self) -> usize {
        return self.n_loops;
    }

    /// the capacity of vertices
    #[inline(always)]
    fn n_vertices_capacity(&self) -> usize {
        return self.v_halfedge_arr.len();
    }

    /// the capacity of halfedges
    #[inline(always)]
    fn n_halfedges_capacity(&self) -> usize {
        return self.he_next_arr.len();
    }

    /// the capacity of edges
    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        return self.e_halfedge_arr.len();
    }

    /// the capacity of faces
    #[inline(always)]
    fn n_faces_capacity(&self) -> usize {
        return self.f_halfedge_arr.len();
    }

    /// the capacity of boundary loops
    #[inline(always)]
    fn n_boundary_loop_capacity(&self) -> usize {
        return self.bl_halfedge_arr.len();
    }

    /// vertex id is valid
    #[inline(always)]
    fn vertex_is_valid(&self, vid: VertexId) -> bool {
        return self.v_halfedge_arr[vid] != INVALID_IND;
    }

    /// halfedge id is valid
    #[inline(always)]
    fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool {
        return self.he_next_arr[hid] != INVALID_IND;
    }

    /// edge id is valid
    #[inline(always)]
    fn edge_is_valid(&self, eid: super::EdgeId) -> bool {
        return self.e_halfedge_arr[eid] != INVALID_IND;
    }

    /// face id is valid
    #[inline(always)]
    fn face_is_valid(&self, fid: FaceId) -> bool {
        return self.f_halfedge_arr[fid] != INVALID_IND;
    }

    /// boundary loop id is valid
    #[inline(always)]
    fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool {
        return self.bl_halfedge_arr[blid] != INVALID_IND;
    }

    /// start vertex iterator
    #[inline(always)]
    fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self> {
        VertexIter::new(vid, self)
    }

    /// start halfedge iterator
    #[inline(always)]
    fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self> {
        HalfedgeIter::new(hid, self)
    }

    /// start edge iterator
    #[inline(always)]
    fn edge<'a>(&'a self, eid: EdgeId) -> super::EdgeIter<'a, Self> {
        EdgeIter::new(eid, self)
    }

    /// start face iterator
    #[inline(always)]
    fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self> {
        FaceIter::new(fid, self)
    }

    /// start boundary loop iterator
    #[inline(always)]
    fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self> {
        BoundaryLoopIter::new(blid, self)
    }

    /// the halfedge starting from this vertex
    #[inline(always)]
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
        HalfedgeId::from(self.v_halfedge_arr[vid])
    }

    /// the start vertex of the halfedge
    #[inline(always)]
    fn he_vertex(&self, hid: HalfedgeId) -> VertexId {
        VertexId::from(self.he_vertex_arr[hid])
    }

    /// the end vertex of the halfedge
    #[inline(always)]
    fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId {
        VertexId::from(self.he_vertex_arr[self.he_next_arr[hid]])
    }

    /// the next halfedge of the halfedge
    #[inline(always)]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_next_arr[hid])
    }

    /// the previous halfedge of the halfedge
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
        let mut curr = hid;
        loop {
            let next = self.he_next(curr);
            if next == hid {
                return curr;
            }
            curr = next;
        }
    }

    /// the twin halfedge of the halfedge
    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_sibling_arr[hid])
    }

    /// the sibling halfedge of the halfedge
    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_sibling_arr[hid])
    }

    /// the next incoming halfedge of the halfedge
    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_vert_in_next_arr[hid])
    }

    /// the next outgoing halfedge of the halfedge
    #[inline(always)]
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_vert_out_next_arr[hid])
    }

    /// the edge of this halfedge
    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        EdgeId::from(self.he_edge_arr[hid])
    }

    /// the face or boundary loop of the halfedge
    #[inline(always)]
    fn he_face_or_boundary_loop(&self, hid: HalfedgeId) -> FaceOrBoundaryLoopId {
        let id = self.he_face_arr[hid];
        if id == INVALID_IND {
            FaceOrBoundaryLoopId::INVALID
        } else if id < self.he_face_arr.len() {
            FaceOrBoundaryLoopId::Face(FaceId::from(id))
        } else {
            FaceOrBoundaryLoopId::BoundaryLoop(BoundaryLoopId::from(BL_START - id))
        }
    }

    /// the first halfedge of the edge
    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        HalfedgeId::from(self.e_halfedge_arr[eid])
    }

    /// the first halfedge of the face
    #[inline(always)]
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
        HalfedgeId::from(self.f_halfedge_arr[fid])
    }

    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        true
    }
}
