use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId, ElementIndex, Halfedge, HalfedgeId};
use crate::{element_id, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

element_id!(struct EdgeId);

pub struct Edge<'a, M: Mesh> {
    mesh: &'a M,
    id: EdgeId,
}

impl<'a, M: Mesh> Edge<'a, M> {
    pub fn new(mesh: &'a M, id: EdgeId) -> Self {
        Self { mesh, id }
    }

    pub fn halfedge(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.e_halfedge(self.id))
    }

    pub fn halfedges(&self) -> EHIter<'a, M> {
        EHIter::new(self.mesh, self.id)
    }
}

impl<'a, M: Mesh> Deref for Edge<'a, M> {
    type Target = EdgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}
pub struct EdgeIter<'a, M: Mesh> {
    mesh: &'a M,
    id: EdgeId,
}

impl<'a, M: Mesh> EdgeIter<'a, M> {
    pub fn new(mesh: &'a M) -> Self {
        let mut iter = Self {
            mesh,
            id: EdgeId(0),
        };
        while !iter.is_end() && !iter.valid() {
            <Self as Element>::next(&mut iter);
        }
        iter
    }
}

impl<'a, M: Mesh> Element for EdgeIter<'a, M> {
    type Item = Edge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Self::Item {
            mesh: self.mesh,
            id: self.id,
        }
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.mesh.edge_is_valid(self.id)
    }

    #[inline(always)]
    fn next(&mut self) {
        *self.id += 1;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        *self.id == self.mesh.n_edges_capacity()
    }
}

impl<'a, M: Mesh> Iterator for EdgeIter<'a, M> {
    type Item = Edge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct EHIter<'a, M: Mesh> {
    mesh: &'a M,
    first: HalfedgeId,
    curr: HalfedgeId,
    just_started: bool,
}

impl<'a, M: Mesh> EHIter<'a, M> {
    pub fn new(mesh: &'a M, eid: EdgeId) -> Self {
        let first = mesh.e_halfedge(eid);
        Self {
            mesh,
            first,
            curr: first,
            just_started: true,
        }
    }
}

impl<'a, M: Mesh> Element for EHIter<'a, M> {
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
        self.just_started = false;
        self.curr = self.mesh.he_sibling(self.curr);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.curr == self.first && !self.just_started
    }
}

impl<'a, M: Mesh> Iterator for EHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
