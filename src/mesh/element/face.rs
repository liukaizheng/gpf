use std::ops::{Deref, Index, IndexMut};

use super::{iter_next, Element};
use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

impl From<usize> for FaceId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl<T> Index<FaceId> for Vec<T> {
    type Output = T;
    fn index(&self, index: FaceId) -> &Self::Output {
        &self[index.0]
    }
}

impl<T> IndexMut<FaceId> for Vec<T> {
    fn index_mut(&mut self, index: FaceId) -> &mut Self::Output {
        &mut self[index.0]
    }
}

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
    type Item = FaceId;
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

    fn next(&mut self) -> bool {
        self.id.0 += 1;
        self.id.0 < self.capacity()
    }
}

impl<'a, M: Mesh> Iterator for FaceIter<'a, M> {
    type Item = FaceId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}