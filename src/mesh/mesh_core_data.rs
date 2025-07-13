use std::alloc::Allocator;

use super::{
    ElementId, FaceId, clone_vec_in,
    element::{HalfedgeId, VertexId},
};

pub struct MeshCoreData<A: Allocator + Copy> {
    pub(crate) v_halfedge_arr: Vec<HalfedgeId, A>,
    pub(crate) he_prev_arr: Vec<HalfedgeId, A>,
    pub(crate) he_next_arr: Vec<HalfedgeId, A>,
    pub(crate) he_vertex_arr: Vec<VertexId, A>,
    pub(crate) he_face_arr: Vec<FaceId, A>,
    pub(crate) f_halfedge_arr: Vec<HalfedgeId, A>,

    pub(crate) n_vertices: usize,
    pub(crate) n_halfedges: usize,
    pub(crate) n_faces: usize,
    pub(crate) alloc: A,
}

impl<A: Allocator + Copy> MeshCoreData<A> {
    pub(crate) fn new(n_vertices: usize, n_faces: usize, alloc: A) -> Self {
        let mut v_halfedge_arr = Vec::with_capacity_in(n_vertices, alloc);
        v_halfedge_arr.resize(n_vertices, HalfedgeId::default());
        let mut f_halfedge_arr = Vec::with_capacity_in(n_faces, alloc);
        f_halfedge_arr.resize(n_faces, HalfedgeId::default());

        Self {
            v_halfedge_arr,
            he_prev_arr: Vec::new_in(alloc),
            he_next_arr: Vec::new_in(alloc),
            he_vertex_arr: Vec::new_in(alloc),
            he_face_arr: Vec::new_in(alloc),
            f_halfedge_arr,
            n_vertices,
            n_halfedges: 0,
            n_faces,
            alloc,
        }
    }

    pub(crate) fn clone_in<A1: Allocator + Copy>(&self, alloc: A1) -> MeshCoreData<A1> {
        MeshCoreData {
            v_halfedge_arr: clone_vec_in(&self.v_halfedge_arr, alloc),
            he_prev_arr: clone_vec_in(&self.he_prev_arr, alloc),
            he_next_arr: clone_vec_in(&self.he_next_arr, alloc),
            he_vertex_arr: clone_vec_in(&self.he_vertex_arr, alloc),
            he_face_arr: clone_vec_in(&self.he_face_arr, alloc),
            f_halfedge_arr: clone_vec_in(&self.f_halfedge_arr, alloc),
            n_vertices: self.n_vertices,
            n_halfedges: self.n_halfedges,
            n_faces: self.n_faces,
            alloc,
        }
    }

    #[inline]
    pub(crate) fn n_vertices_capacity(&self) -> usize {
        self.v_halfedge_arr.len()
    }

    #[inline]
    pub(crate) fn n_halfedges_capacity(&self) -> usize {
        self.he_vertex_arr.len()
    }

    #[inline]
    pub(crate) fn n_faces_capacity(&self) -> usize {
        self.f_halfedge_arr.len()
    }

    #[inline]
    pub(crate) fn connect_halfedges(&mut self, hid0: HalfedgeId, hid1: HalfedgeId) {
        self.he_next_arr[hid0] = hid1;
        self.he_prev_arr[hid1] = hid0;
    }

    #[inline]
    pub(crate) fn set_v_halfedge(&mut self, vid: VertexId, hid: HalfedgeId) {
        if vid.valid() {
            self.v_halfedge_arr[vid] = hid;
        }
    }

    #[inline]
    pub(crate) fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId) {
        self.he_vertex_arr[hid] = vid;
    }

    #[inline]
    pub(crate) fn set_he_face(&mut self, hid: HalfedgeId, fid: FaceId) {
        self.he_face_arr[hid] = fid;
    }

    #[inline]
    pub(crate) fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId) {
        self.f_halfedge_arr[fid] = hid;
    }

    #[inline]
    pub(crate) fn v_min_reserve(&mut self, vid: VertexId) {
        let len = vid.0 + 1;
        if self.v_halfedge_arr.len() < len {
            self.v_halfedge_arr.resize(len, HalfedgeId::default());
        }
    }

    #[inline]
    pub(crate) fn recount_n_vertices(&mut self) {
        self.n_vertices = self.v_halfedge_arr.len();
        for &hid in &self.v_halfedge_arr {
            if !hid.valid() {
                self.n_vertices -= 1;
            }
        }
    }

    #[inline]
    pub(crate) fn reserve_halfedges(&mut self, n_halfedges: usize) {
        self.he_prev_arr.reserve(n_halfedges);
        self.he_next_arr.reserve(n_halfedges);
        self.he_vertex_arr.reserve(n_halfedges);
        self.he_face_arr.reserve(n_halfedges);
    }

    #[inline]
    pub(crate) fn new_halfedges(&mut self, n_halfedges: usize) {
        let new_len = self.he_vertex_arr.len() + n_halfedges;
        self.he_prev_arr.resize(new_len, HalfedgeId::default());
        self.he_next_arr.resize(new_len, HalfedgeId::default());
        self.he_vertex_arr.resize(new_len, VertexId::default());
        self.he_face_arr.resize(new_len, FaceId::default());
    }
}
