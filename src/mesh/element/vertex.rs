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
        VVIter {
            imp: VNeiIterImpl::new(self.mesh, self.id),
        }
    }

    pub fn edges(&self) -> VEIter<'a, M> {
        VEIter {
            imp: VNeiIterImpl::new(self.mesh, self.id),
        }
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
    he: HalfedgeId,
    incoming: bool,
}

impl VNeiState {
    fn new(he: HalfedgeId, incoming: bool) -> Self {
        Self { he, incoming }
    }

    #[inline]
    fn is_halfedge_canonical<M: Mesh>(&self, mesh: &M) -> bool {
        if mesh.use_implicit_twin() {
            true
        } else {
            let edge = mesh.he_edge(self.he);
            self.he == mesh.e_halfedge(edge)
        }
    }

    fn next<M: Mesh>(&mut self, mesh: &M, first_he: &mut HalfedgeId) {
        if mesh.use_implicit_twin() {
            self.he = mesh.he_next_outgoing_neighbor(self.he);
        } else {
            if self.incoming {
                self.he = mesh.he_next_incoming_neighbor(self.he);
                if self.he == *first_he {
                    self.incoming = false;
                    self.he = mesh.he_next(self.he);
                    *first_he = self.he;
                }
            } else {
                self.he = mesh.he_next_outgoing_neighbor(self.he);
                if self.he == *first_he {
                    self.incoming = true;
                    self.he = mesh.he_prev(self.he);
                    *first_he = self.he;
                }
            }
        }
    }
}

struct VNeiIterImpl<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    curr: VNeiState,
    start: VNeiState,
    first_he: HalfedgeId,
}

impl<'a, M: Mesh> VNeiIterImpl<'a, M> {
    fn new(mesh: &'a M, id: VertexId) -> Self {
        let mut first_he = mesh.v_halfedge(id);
        let mut curr = VNeiState::new(first_he, false);
        while !curr.is_halfedge_canonical(mesh) {
            curr.next(mesh, &mut first_he);
        }
        let start = curr.clone();
        first_he = start.he;
        Self {
            mesh,
            just_start: true,
            curr,
            start,
            first_he,
        }
    }
}

pub struct VVIter<'a, M: Mesh> {
    imp: VNeiIterImpl<'a, M>,
}

impl<'a, M: Mesh> Element for VVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        let vid = if self.imp.curr.incoming {
            self.imp.mesh.he_vertex(self.imp.curr.he)
        } else {
            self.imp.mesh.he_tip_vertex(self.imp.curr.he)
        };
        Vertex::new(self.imp.mesh, vid)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.imp.curr.is_halfedge_canonical(self.imp.mesh)
    }

    #[inline(always)]
    fn next(&mut self) {
        self.imp.just_start = false;
        self.imp.curr.next(self.imp.mesh, &mut self.imp.first_he);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.imp.curr == self.imp.start && !self.imp.just_start
    }
}

impl<'a, M: Mesh> Iterator for VVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct VEIter<'a, M: Mesh> {
    imp: VNeiIterImpl<'a, M>,
}

impl<'a, M: Mesh> Element for VEIter<'a, M> {
    type Item = Edge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        let eid = self.imp.mesh.he_edge(self.imp.curr.he);
        Edge::new(self.imp.mesh, eid)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.imp.curr.is_halfedge_canonical(self.imp.mesh)
    }

    #[inline(always)]
    fn next(&mut self) {
        self.imp.just_start = false;
        self.imp.curr.next(self.imp.mesh, &mut self.imp.first_he);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.imp.curr == self.imp.start && !self.imp.just_start
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
        let start = mesh.v_halfedge(id);
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
        Halfedge::new(self.mesh, self.curr)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.just_start = false;
        self.curr = self.mesh.he_next_outgoing_neighbor(self.curr);
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
