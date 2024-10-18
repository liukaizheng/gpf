use super::element::{HalfedgeId, VertexId};

pub struct MeshCoreData {
    pub(crate) v_halfedge_arr: Vec<HalfedgeId>,
    pub(crate) he_prev_arr: Vec<HalfedgeId>,
    pub(crate) he_next_arr: Vec<HalfedgeId>,
    pub(crate) he_vertex_arr: Vec<VertexId>,
    pub(crate) f_halfedge_arr: Vec<HalfedgeId>,

    pub(crate) n_vertices: usize,
    pub(crate) n_halfedges: usize,
    pub(crate) n_faces: usize,
}

impl MeshCoreData {
    pub(crate) fn new(n_vertices: usize, n_faces: usize) -> Self {
        Self {
            v_halfedge_arr: vec![HalfedgeId::default(); n_vertices],
            he_prev_arr: Vec::new(),
            he_next_arr: Vec::new(),
            he_vertex_arr: Vec::new(),
            f_halfedge_arr: vec![HalfedgeId::default(); n_faces],
            n_vertices,
            n_halfedges: 0,
            n_faces,
        }
    }

    #[inline]
    pub(crate) fn connect_halfedges(&mut self, hid0: HalfedgeId, hid1: HalfedgeId) {
        self.he_next_arr[hid0] = hid1;
        self.he_prev_arr[hid1] = hid0;
    }
}
