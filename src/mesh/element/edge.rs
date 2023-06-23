use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId, HalfedgeId, HalfedgeIter};
use crate::{element_id, element_iterator, halfedges_iterator, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

element_id! {struct EdgeId}

pub struct EdgeIter<'a, M: Mesh> {
    id: EdgeId,
    mesh: &'a M,
}

impl<'a, M: Mesh> EdgeIter<'a, M> {
    pub fn new(id: EdgeId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_edges()
    }

    #[inline(always)]
    fn capacity(&self) -> usize {
        self.mesh.n_edges_capacity()
    }
}

impl<'a, M: Mesh> Deref for EdgeIter<'a, M> {
    type Target = EdgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

element_iterator! {
    struct EdgeIter -> EdgeId, {
        fn id(&self) -> EdgeId {
            self.id
        }
    }, {
        fn valid(&self) -> bool {
            self.mesh.edge_is_valid(self.id)
        }
    },{
        fn next(&mut self) {
            *self.id += 1;
        }
    },{
        fn is_end(&self) -> bool {
            *self.id == self.capacity()
        }
    }
}

pub trait Edge: Element {
    fn halfedge(&self) -> HalfedgeIter<Self::M>;
    fn halfedges(&self) -> EHalfedgesIter<Self::M>;
}

impl<'a, M: Mesh> Edge for EdgeIter<'a, M> {
    fn halfedge(&self) -> HalfedgeIter<Self::M> {
        HalfedgeIter::new(self.mesh.e_halfedge(self.id), self.mesh)
    }
    fn halfedges(&self) -> EHalfedgesIter<Self::M> {
        EHalfedgesIter::new(self.mesh.e_halfedge(self.id), self.mesh)
    }
}

pub struct EHalfedgesIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, M: Mesh> EHalfedgesIter<'a, M> {
    pub fn new(he: HalfedgeId, mesh: &'a M) -> Self {
        Self {
            mesh,
            just_start: true,
            first_he: he,
            curr_he: he,
        }
    }

    fn next(&mut self) {
        self.curr_he = self.mesh.he_sibling(self.curr_he);
    }
}

halfedges_iterator! {struct EHalfedgesIter}
