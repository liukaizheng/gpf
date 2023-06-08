use std::{
    marker::PhantomData,
    ops::{Deref, DerefMut, Index, IndexMut},
};

use super::{iter_next, Element, ElementId, HalfedgeId, HalfedgeIter};
use crate::{element_id, element_iterator, halfedges_iterator, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

element_id! {struct EdgeId}

pub struct EdgeIter<'a, 'b: 'a, M: Mesh<'b>> {
    phantom: PhantomData<&'b M>,
    id: EdgeId,
    mesh: &'a M,
}

impl<'a, 'b: 'a, M: Mesh<'b>> EdgeIter<'a, 'b, M> {
    pub fn new(id: EdgeId, mesh: &'a M) -> Self {
        Self {
            id,
            mesh,
            phantom: PhantomData,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_halfedges()
    }

    #[inline(always)]
    fn capacity(&self) -> usize {
        self.mesh.n_halfedges_capacity()
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Deref for EdgeIter<'a, 'b, M> {
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

pub trait Edge<'b>: Element<'b> {
    fn halfedge(&self) -> HalfedgeIter<'_, 'b, Self::M>;
    fn halfedges(&self) -> EHalfedgesIter<'_, 'b, Self::M>;
}

impl<'a, 'b: 'a, M: Mesh<'b>> Edge<'b> for EdgeIter<'a, 'b, M> {
    fn halfedge(&self) -> HalfedgeIter<'_, 'b, Self::M> {
        HalfedgeIter::new(self.mesh.e_halfedge(self.id), self.mesh)
    }
    fn halfedges(&self) -> EHalfedgesIter<'_, 'b, Self::M> {
        EHalfedgesIter::new(self.mesh.e_halfedge(self.id), self.mesh)
    }
}

pub struct EHalfedgesIter<'a, 'b: 'a, M: Mesh<'b>> {
    phantom: PhantomData<&'b M>,
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, 'b: 'a, M: Mesh<'b>> EHalfedgesIter<'a, 'b, M> {
    pub fn new(he: HalfedgeId, mesh: &'a M) -> Self {
        Self {
            phantom: PhantomData,
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
