use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, EdgeIter, Element, ElementId};
use crate::{element_id, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(pub usize);

element_id! {struct HalfedgeId}

pub struct HalfedgeIter<'a, M: Mesh> {
    id: HalfedgeId,
    mesh: &'a M,
}

impl<'a, M: Mesh> HalfedgeIter<'a, M> {
    pub fn new(id: HalfedgeId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_halfedges()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_halfedges_capacity()
    }
}

impl<'a, M: Mesh> Deref for HalfedgeIter<'a, M> {
    type Target = HalfedgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, M: Mesh> Element for HalfedgeIter<'a, M> {
    type Item = HalfedgeId;
    type M = M;

    fn id(&self) -> HalfedgeId {
        self.id
    }

    fn mesh(&self) -> &M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.halfedge_is_valid(self.id)
    }

    fn next(&mut self) {
        *self.id += 1;
    }

    fn is_end(&self) -> bool {
        *self.id == self.capacity()
    }
}

pub trait Halfedge: Element {
    fn edge(&self) -> EdgeIter<'_, Self::M>;
    fn next(&self) -> HalfedgeIter<'_, Self::M>;
    fn prev(&self) -> HalfedgeIter<'_, Self::M>;
}

impl<'a, M: Mesh> Halfedge for HalfedgeIter<'a, M> {
    fn edge(&self) -> EdgeIter<'_, Self::M> {
        EdgeIter::new(self.mesh.he_edge(self.id), self.mesh)
    }
    fn next(&self) -> HalfedgeIter<'_, Self::M> {
        HalfedgeIter::new(self.mesh.he_next(self.id), self.mesh)
    }
    fn prev(&self) -> HalfedgeIter<'_, Self::M> {
        HalfedgeIter::new(self.mesh.he_prev(self.id), self.mesh)
    }
}

impl<'a, M: Mesh> Iterator for HalfedgeIter<'a, M> {
    type Item = HalfedgeId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
