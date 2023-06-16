use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId, HalfedgeId, HalfedgeIter};
use crate::{element_id, element_iterator, halfedges_iterator, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

element_id! {struct FaceId}

pub struct FaceIter<'a, M: Mesh> {
    id: FaceId,
    mesh: &'a M,
}

impl<'a, M: Mesh> FaceIter<'a, M> {
    pub fn new(id: FaceId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_faces()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_faces_capacity()
    }
}

impl<'a, M: Mesh> Deref for FaceIter<'a, M> {
    type Target = FaceId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, M: Mesh> Element for FaceIter<'a, M> {
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

impl<'a, M: Mesh> Iterator for FaceIter<'a, M> {
    type Item = FaceId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
pub trait Face: Element {
    fn halfedge(&self) -> HalfedgeIter<Self::M>;
    fn halfedges(&self) -> FHalfedgesIter<Self::M>;
}

impl<'a, M: Mesh> Face for FaceIter<'a, M> {
    fn halfedge(&self) -> HalfedgeIter<Self::M> {
        HalfedgeIter::new(self.mesh.f_halfedge(self.id), self.mesh)
    }
    fn halfedges(&self) -> FHalfedgesIter<Self::M> {
        FHalfedgesIter::new(self.mesh.f_halfedge(self.id), self.mesh)
    }
}

pub struct FHalfedgesIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, M: Mesh> FHalfedgesIter<'a, M> {
    pub fn new(he: HalfedgeId, mesh: &'a M) -> Self {
        Self {
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
