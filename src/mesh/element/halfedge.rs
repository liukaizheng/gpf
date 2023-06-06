use std::{ops::{Deref, DerefMut, Index, IndexMut}, marker::PhantomData};

use super::{iter_next, EdgeIter, Element, ElementId};
use crate::{element_id, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(pub usize);

element_id! {struct HalfedgeId}

pub struct HalfedgeIter<'a, 'b: 'a, M: Mesh<'b>> {
    phantom: PhantomData<&'b M>,
    id: HalfedgeId,
    mesh: &'a M,
}

impl<'a, 'b: 'a, M: Mesh<'b>> HalfedgeIter<'a, 'b, M> {
    pub fn new(id: HalfedgeId, mesh: &'a M) -> Self {
        Self { id, mesh, phantom: PhantomData}
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_halfedges()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_halfedges_capacity()
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Deref for HalfedgeIter<'a, 'b, M> {
    type Target = HalfedgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Element<'b> for HalfedgeIter<'a, 'b, M> {
    type Id = HalfedgeId;
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

pub trait Halfedge<'b>: Element<'b> {
    fn edge(&self) -> EdgeIter<'_, 'b, Self::M>;
    fn next(&self) -> HalfedgeIter<'_, 'b, Self::M>;
    fn prev(&self) -> HalfedgeIter<'_, 'b, Self::M>;
}

impl<'a, 'b: 'a, M: Mesh<'b>> Halfedge<'b> for HalfedgeIter<'a, 'b, M> {
    fn edge(&self) -> EdgeIter<'_, 'b, Self::M> {
        EdgeIter::new(self.mesh.he_edge(self.id), self.mesh)
    }
    fn next(&self) -> HalfedgeIter<'_, 'b, Self::M> {
        HalfedgeIter::new(self.mesh.he_next(self.id), self.mesh)
    }
    fn prev(&self) -> HalfedgeIter<'_, 'b, Self::M> {
        HalfedgeIter::new(self.mesh.he_prev(self.id), self.mesh)
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Iterator for HalfedgeIter<'a, 'b, M> {
    type Item = HalfedgeId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
