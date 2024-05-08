use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId, ElementIndex, Halfedge, HalfedgeId, Vertex};
use crate::{element_id, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

element_id!(struct FaceId);

pub struct Face<'a, M: Mesh> {
    mesh: &'a M,
    id: FaceId,
}

impl<'a, M: Mesh> Face<'a, M> {
    pub fn new(mesh: &'a M, id: FaceId) -> Self {
        Self { mesh, id }
    }

    pub fn halfedge(&self) -> Halfedge<'a, M> {
        Halfedge::new(self.mesh, self.mesh.f_halfedge(self.id))
    }

    pub fn halfedges(&self) -> FHIter<'a, M> {
        FHIter::new(self.mesh, self.id)
    }

    pub fn vertices(&self) -> FVIter<'a, M> {
        FVIter::new(self.mesh, self.id)
    }
}

impl<'a, M: Mesh> Deref for Face<'a, M> {
    type Target = FaceId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

pub struct FaceIter<'a, M: Mesh> {
    mesh: &'a M,
    id: FaceId,
}

impl<'a, M: Mesh> FaceIter<'a, M> {
    pub fn new(mesh: &'a M) -> Self {
        let mut iter = Self {
            mesh,
            id: FaceId(0),
        };
        while !iter.is_end() && !iter.valid() {
            <Self as Element>::next(&mut iter);
        }
        iter
    }
}

impl<'a, M: Mesh> Element for FaceIter<'a, M> {
    type Item = Face<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Self::Item {
            mesh: self.mesh,
            id: self.id,
        }
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        self.mesh.face_is_valid(self.id)
    }

    #[inline(always)]
    fn next(&mut self) {
        *self.id += 1;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        *self.id == self.mesh.n_faces_capacity()
    }
}

impl<'a, M: Mesh> Iterator for FaceIter<'a, M> {
    type Item = Face<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct FHIter<'a, M: Mesh> {
    mesh: &'a M,
    first: HalfedgeId,
    current: HalfedgeId,
    just_start: bool,
}

impl<'a, M: Mesh> FHIter<'a, M> {
    pub fn new(mesh: &'a M, fid: FaceId) -> Self {
        let first = mesh.f_halfedge(fid);
        Self {
            mesh,
            first,
            current: first,
            just_start: true,
        }
    }
}

impl<'a, M: Mesh> Element for FHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.mesh, self.current)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.current = self.mesh.he_next(self.current);
        self.just_start = false;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.current == self.first && !self.just_start
    }
}

impl<'a, M: Mesh> Iterator for FHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

pub struct FVIter<'a, M: Mesh> {
    mesh: &'a M,
    he: HalfedgeId,
    start: HalfedgeId,
    just_start: bool,
}

impl<'a, M: Mesh> FVIter<'a, M> {
    pub fn new(mesh: &'a M, fid: FaceId) -> Self {
        let he = mesh.f_halfedge(fid);
        Self {
            mesh,
            he,
            start: he,
            just_start: true,
        }
    }
}

impl<'a, M: Mesh> Element for FVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Vertex::new(self.mesh, self.mesh.he_vertex(self.he))
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.he = self.mesh.he_next(self.he);
        self.just_start = false;
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.he == self.start && !self.just_start
    }
}

impl<'a, M: Mesh> Iterator for FVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
