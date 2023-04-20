use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId};
use crate::{element_id, mesh::Mesh, INVALID_IND};

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
