use std::{
    marker::PhantomData,
    ops::{Deref, DerefMut, Index, IndexMut},
};

use super::{iter_next, Element, ElementId, HalfedgeId, HalfedgeIter};
use crate::{element_id, element_iterator, halfedges_iterator, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

element_id! {struct FaceId}

pub struct FaceIter<'a, 'b: 'a, M: Mesh<'b>> {
    phantom: PhantomData<&'b M>,
    id: FaceId,
    mesh: &'a M,
}

impl<'a, 'b: 'a, M: Mesh<'b>> FaceIter<'a, 'b, M> {
    pub fn new(id: FaceId, mesh: &'a M) -> Self {
        Self {
            id,
            mesh,
            phantom: PhantomData,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_faces()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_faces_capacity()
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Deref for FaceIter<'a, 'b, M> {
    type Target = FaceId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Element<'b> for FaceIter<'a, 'b, M> {
    type Id = FaceId;
    type M = M;

    fn id(&self) -> FaceId {
        self.id
    }

    fn mesh(&self) -> &M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.face_is_valid(self.id)
    }

    fn next(&mut self) {
        self.id.0 += 1;
    }

    fn is_end(&self) -> bool {
        *self.id == self.capacity()
    }
}

impl<'a, 'b: 'a, M: Mesh<'b>> Iterator for FaceIter<'a, 'b, M> {
    type Item = FaceId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
pub trait Face<'b>: Element<'b> {
    fn halfedge(&self) -> HalfedgeIter<'_, 'b, Self::M>;
    fn halfedges(&self) -> FHalfedgesIter<'_, 'b, Self::M>;
}

impl<'a, 'b: 'a, M: Mesh<'b>> Face<'b> for FaceIter<'a, 'b, M> {
    fn halfedge(&self) -> HalfedgeIter<'_, 'b, Self::M> {
        HalfedgeIter::new(self.mesh.f_halfedge(self.id), self.mesh)
    }
    fn halfedges(&self) -> FHalfedgesIter<'_, 'b, Self::M> {
        FHalfedgesIter::new(self.mesh.f_halfedge(self.id), self.mesh)
    }
}

pub struct FHalfedgesIter<'a, 'b: 'a, M: Mesh<'b>> {
    phantom: PhantomData<&'b M>,
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, 'b: 'a, M: Mesh<'b>> FHalfedgesIter<'a, 'b, M> {
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
        self.curr_he = self.mesh.he_next(self.curr_he);
    }
}

halfedges_iterator! {struct FHalfedgesIter}
