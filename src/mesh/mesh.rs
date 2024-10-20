use super::element::{
    Edge, EdgeId, EdgeIter, ElementId, ElementIndex, Face, FaceId, FaceIter, Halfedge, HalfedgeId,
    HalfedgeIter, Vertex, VertexId, VertexIter,
};
use super::mesh_core_data::MeshCoreData;

pub trait Mesh: Sized {
    fn core_data(&self) -> &MeshCoreData;

    #[inline(always)]
    fn n_vertices(&self) -> usize {
        self.core_data().n_vertices
    }

    #[inline(always)]
    fn n_halfedges(&self) -> usize {
        self.core_data().n_halfedges
    }
    fn n_edges(&self) -> usize;

    fn n_faces(&self) -> usize {
        self.core_data().n_faces
    }
    // fn n_boundary_loops(&self) -> usize;

    #[inline(always)]
    fn n_vertices_capacity(&self) -> usize {
        self.core_data().v_halfedge_arr.len()
    }

    #[inline(always)]
    fn n_halfedges_capacity(&self) -> usize {
        self.core_data().he_next_arr.len()
    }

    fn n_edges_capacity(&self) -> usize;

    #[inline(always)]
    fn n_faces_capacity(&self) -> usize {
        self.core_data().f_halfedge_arr.len()
    }

    #[inline(always)]
    fn vertex_is_valid(&self, vid: VertexId) -> bool {
        self.v_halfedge(vid).valid()
    }

    #[inline(always)]
    fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool {
        self.he_vertex(hid).valid()
    }

    fn edge_is_valid(&self, eid: EdgeId) -> bool;

    #[inline(always)]
    fn face_is_valid(&self, fid: FaceId) -> bool {
        self.f_halfedge(fid).valid()
    }

    #[inline(always)]
    fn vertex(&self, vid: VertexId) -> Vertex<Self> {
        Vertex::new(self, vid)
    }

    #[inline(always)]
    fn halfedge(&self, hid: HalfedgeId) -> Halfedge<Self> {
        Halfedge::new(self, hid)
    }

    #[inline(always)]
    fn edge(&self, eid: EdgeId) -> Edge<Self> {
        Edge::new(self, eid)
    }

    #[inline(always)]
    fn face(&self, fid: FaceId) -> Face<Self> {
        Face::new(self, fid)
    }

    #[inline(always)]
    fn vertices(&self) -> VertexIter<Self> {
        VertexIter::new(self)
    }

    #[inline(always)]
    fn halfedges(&self) -> HalfedgeIter<Self> {
        HalfedgeIter::new(self)
    }

    #[inline(always)]
    fn edges(&self) -> EdgeIter<Self> {
        EdgeIter::new(self)
    }

    #[inline(always)]
    fn faces(&self) -> FaceIter<Self> {
        FaceIter::new(self)
    }

    /// the halfedge starting from this vertex
    #[inline(always)]
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
        self.core_data().v_halfedge_arr[vid.index()]
    }

    /// the start vertex of the halfedge
    #[inline(always)]
    fn he_vertex(&self, hid: HalfedgeId) -> VertexId {
        self.core_data().he_vertex_arr[hid.index()]
    }
    /// the end vertex of the halfedge
    fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId {
        self.he_vertex(self.he_next(hid))
    }

    /// the next halfedge of the halfedge
    #[inline(always)]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.core_data().he_next_arr[hid.index()]
    }

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

    /// the first halfedge of the edge
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId;

    /// two vertices of the edge
    #[inline(always)]
    fn e_vertices(&self, eid: EdgeId) -> [VertexId; 2] {
        let hid = self.e_halfedge(eid);
        [self.he_vertex(hid), self.he_tip_vertex(hid)]
    }

    /// the first halfedge of the face
    #[inline(always)]
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
        self.core_data().f_halfedge_arr[fid.index()]
    }

    #[inline(always)]
    fn he_same_dir(&self, hid: HalfedgeId) -> bool {
        let eid = self.he_edge(hid);
        self.he_vertex(hid) == self.he_vertex(self.e_halfedge(eid))
    }

    #[inline(always)]
    fn hes_same_dir(&self, ha: HalfedgeId, hb: HalfedgeId) -> bool {
        self.he_vertex(ha) == self.he_vertex(hb)
    }

    fn use_implicit_twin(&self) -> bool;
}
