use std::{cell::RefCell, rc::Weak};

use super::{
    mesh_core_data::MeshCoreData, BoundaryLoopId, EdgeId, HalfedgeId, Mesh, VertexId
};

pub struct ManifoldMesh {
    core_data: MeshCoreData,

    he_face_arr: Vec<usize>,
    v_halfedge_arr: Vec<HalfedgeId>,
    he_next_arr: Vec<HalfedgeId>,
    he_vertex_arr: Vec<VertexId>,

}

impl ManifoldMesh {


}

impl Mesh for ManifoldMesh {
    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        true
    }

    #[inline(always)]
    fn core_data(&self) -> &MeshCoreData {
        &self.core_data
    }

    #[inline(always)]
    fn n_edges(&self) -> usize {
        self.core_data.n_halfedges >> 1
    }

    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        self.n_halfedges_capacity() >> 1
    }

    #[inline(always)]
    fn edge_is_valid(&self, eid: EdgeId) -> bool {
        let idx = eid.0 << 1;
        self.halfedge_is_valid(idx.into()) || self.halfedge_is_valid((idx + 1).into())
    }

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

    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        (hid.0 ^ 1).into()
    }

    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_twin(hid)
    }

    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_vert_in_next_arr[hid]
    }

    #[inline(always)]
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_vert_out_next_arr[hid]
    }

    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        (hid.0 >> 1).into()
    }

    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        let idx = eid.0 << 1;
        if self.halfedge_is_valid(idx.into()) {
            idx.into()
        } else {
            (idx + 1).into()
        }
    }
}