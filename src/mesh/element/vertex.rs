use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Edge, Element, ElementId, ElementIndex, Halfedge, HalfedgeId, INVALID_IND};
use crate::element_id;
use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(pub usize);

element_id!(struct VertexId);
pub struct Vertex<'a, M: Mesh> {
    mesh: &'a M,
    id: VertexId,
}

impl<'a, M: Mesh> Vertex<'a, M> {
    pub fn new(mesh: &'a M, id: VertexId) -> Self {
        Self { mesh, id }
    }
    pub fn in_halfedge(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.he_prev(self.mesh.v_halfedge(self.id)))
    }

    pub fn out_halfedge(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.v_halfedge(self.id))
    }

    pub fn vertices(&self) -> VVIter<'a, M> {
        VVIter::new(self.mesh, self.id)
    }

    pub fn edges(&self) -> VEIter<'a, M> {
        VEIter::new(self.mesh, self.id)
    }

    pub fn incoming_halfedges(&self) -> VIHIter<'a, M> {
        VIHIter::new(self.mesh, self.id)
    }

    pub fn outgoing_halfedges(&self) -> VOHIter<'a, M> {
        VOHIter::new(self.mesh, self.id)
    }
}

impl<'a, M: Mesh> Deref for Vertex<'a, M> {
    type Target = VertexId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}
pub struct VertexIter<'a, M: Mesh> {
    mesh: &'a M,
    id: VertexId,
}

impl<'a, M: Mesh> VertexIter<'a, M> {
    pub fn new(mesh: &'a M) -> Self {
        let mut iter = Self {
            mesh,
            id: VertexId(0),
        };
        while !iter.is_end() && !iter.valid() {
            <Self as Element>::next(&mut iter);
        }
        iter
    }
}

impl<'a, M: Mesh> Element for VertexIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Self::Item {
            mesh: self.mesh,
            id: self.id,
        }
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.mesh.vertex_is_valid(self.id)
    }

    #[inline(always)]
    fn next(&mut self) {
        *self.id += 1;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        *self.id == self.mesh.n_vertices_capacity()
    }
}

impl<'a, M: Mesh> Iterator for VertexIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

#[derive(Clone, PartialEq, Eq)]
struct VNeiState {
    hid: HalfedgeId,
    prev_visited: bool,
}

impl VNeiState {
    #[inline(always)]
    fn new(hid: HalfedgeId) -> Self {
        Self {
            hid,
            prev_visited: false,
        }
    }

    fn halfedge<M: Mesh>(&self, mesh: &M) -> HalfedgeId {
        if !self.prev_visited {
            self.hid
        } else {
            mesh.he_next(self.hid)
        }
    }

    #[inline(always)]
    fn is_canonical<M: Mesh>(&self, mesh: &M) -> bool {
        let hid = self.halfedge(mesh);
        mesh.e_halfedge(mesh.he_edge(hid)) == hid
    }

    fn next<M: Mesh>(&mut self, mesh: &M) {
        if self.prev_visited {
            self.hid = mesh.he_next_incoming_neighbor(self.hid);
            self.prev_visited = false;
        } else {
            self.prev_visited = true;
        }
    }
}

#[derive(Clone, PartialEq, Eq)]
struct VNeighbor<'a, M: Mesh> {
    mesh: &'a M,
    state: VNeiState,
    first_state: VNeiState,
    first: bool,
}

impl<'a, M: Mesh> VNeighbor<'a, M> {
    fn new(mesh: &'a M, vid: VertexId) -> Self {
        let hid = mesh.he_prev(mesh.v_halfedge(vid));
        let mut state = VNeiState::new(hid);
        loop {
            if state.is_canonical(mesh) {
                break;
            }
            state.next(mesh)
        }
        let first_state = state.clone();
        Self {
            mesh,
            state,
            first_state,
            first: true,
        }
    }

    #[inline(always)]
    fn next(&mut self) {
        self.state.next(self.mesh);
        if self.first {
            self.first = false;
        }
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        !self.first && self.state == self.first_state
    }

    #[inline(always)]
    fn halfedge(&self) -> HalfedgeId {
        self.state.halfedge(self.mesh)
    }

    #[inline(always)]
    fn is_canonical(&self) -> bool {
        self.state.is_canonical(self.mesh)
    }
}

pub struct VVIter<'a, M: Mesh> {
    imp: VNeighbor<'a, M>,
}

impl<'a, M: Mesh> VVIter<'a, M> {
    fn new(mesh: &'a M, vid: VertexId) -> Self {
        Self {
            imp: VNeighbor::new(mesh, vid),
        }
    }
}

impl<'a, M: Mesh> Element for VVIter<'a, M> {
    type Item = Vertex<'a, M>;

    fn item(&self) -> Self::Item {
        let hid = self.imp.halfedge();
        let mesh = self.imp.mesh;
        if self.imp.state.prev_visited {
            Vertex::new(mesh, mesh.he_to(hid))
        } else {
            Vertex::new(mesh, mesh.he_from(hid))
        }
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.imp.is_canonical()
    }

    #[inline(always)]
    fn next(&mut self) {
        self.imp.next()
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.imp.is_end()
    }
}

impl<'a, M: Mesh> Iterator for VVIter<'a, M> {
    type Item = Vertex<'a, M>;

    //#[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct VEIter<'a, M: Mesh> {
    imp: VNeighbor<'a, M>,
}

impl<'a, M: Mesh> VEIter<'a, M> {
    fn new(mesh: &'a M, vid: VertexId) -> Self {
        Self {
            imp: VNeighbor::new(mesh, vid),
        }
    }
}

impl<'a, M: Mesh> Element for VEIter<'a, M> {
    type Item = Edge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        let hid = self.imp.halfedge();
        let mesh = self.imp.mesh;
        Edge::new(mesh, mesh.he_edge(hid))
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.imp.is_canonical()
    }

    #[inline(always)]
    fn next(&mut self) {
        self.imp.next()
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.imp.is_end()
    }
}

impl<'a, M: Mesh> Iterator for VEIter<'a, M> {
    type Item = Edge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct VIHIter<'a, M: Mesh> {
    mesh: &'a M,
    curr: HalfedgeId,
    start: HalfedgeId,
    just_start: bool,
}

impl<'a, M: Mesh> VIHIter<'a, M> {
    pub fn new(mesh: &'a M, vid: VertexId) -> Self {
        let start = mesh.he_prev(mesh.v_halfedge(vid));
        Self {
            mesh,
            curr: start,
            start,
            just_start: true,
        }
    }
}

impl<'a, M: Mesh> Element for VIHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.mesh, self.curr)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.just_start = false;
        self.curr = self.mesh.he_next_incoming_neighbor(self.curr);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.curr == self.start && !self.just_start
    }
}

impl<'a, M: Mesh> Iterator for VIHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct VOHIter<'a, M: Mesh> {
    mesh: &'a M,
    curr: HalfedgeId,
    start: HalfedgeId,
    just_start: bool,
}

impl<'a, M: Mesh> VOHIter<'a, M> {
    pub fn new(mesh: &'a M, id: VertexId) -> Self {
        let start = mesh.he_prev(mesh.v_halfedge(id));
        Self {
            mesh,
            curr: start,
            start,
            just_start: true,
        }
    }
}

impl<'a, M: Mesh> Element for VOHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.mesh, self.mesh.he_next(self.curr))
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.just_start = false;
        self.curr = self.mesh.he_next_incoming_neighbor(self.curr);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.curr == self.start && !self.just_start
    }
}

impl<'a, M: Mesh> Iterator for VOHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
