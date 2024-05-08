use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Edge, Element, ElementId, ElementIndex, Vertex};
use crate::{element_id, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(pub usize);

element_id!(struct HalfedgeId);

pub struct Halfedge<'a, M: Mesh> {
    mesh: &'a M,
    id: HalfedgeId,
}

impl<'a, M: Mesh> Halfedge<'a, M> {
    pub fn new(mesh: &'a M, id: HalfedgeId) -> Self {
        Self { mesh, id }
    }

    pub fn from(&self) -> Vertex<'a, M> {
        Vertex::new(self.mesh, self.mesh.he_vertex(self.id))
    }

    pub fn to(&self) -> Vertex<'a, M> {
        Vertex::new(self.mesh, self.mesh.he_tip_vertex(self.id))
    }

    pub fn vertices(&self) -> [Vertex<'a, M>; 2] {
        [self.from(), self.to()]
    }

    pub fn prev(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.he_prev(self.id))
    }

    pub fn next(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.he_next(self.id))
    }

    pub fn twin(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.he_twin(self.id))
    }

    pub fn sibling(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.he_sibling(self.id))
    }

    pub fn edge(&self) -> Edge<'a, M> {
        Edge::new(self.mesh, self.mesh.he_edge(self.id))
    }
}

impl<'a, M: Mesh> Deref for Halfedge<'a, M> {
    type Target = HalfedgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

pub struct HalfedgeIter<'a, M: Mesh> {
    mesh: &'a M,
    id: HalfedgeId,
}

impl<'a, M: Mesh> HalfedgeIter<'a, M> {
    pub fn new(mesh: &'a M) -> Self {
        let mut iter = Self {
            mesh,
            id: HalfedgeId(0),
        };
        while !iter.is_end() && !iter.valid() {
            <Self as Element>::next(&mut iter);
        }
        iter
    }
}

impl<'a, M: Mesh> Element for HalfedgeIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Self::Item {
            mesh: self.mesh,
            id: self.id,
        }
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.mesh.halfedge_is_valid(self.id)
    }

    #[inline(always)]
    fn next(&mut self) {
        *self.id += 1;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        *self.id == self.mesh.n_halfedges_capacity()
    }
}

impl<'a, M: Mesh> Iterator for HalfedgeIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
