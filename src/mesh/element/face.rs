use std::{ops::{Deref, DerefMut, Index, IndexMut}, marker::PhantomData};

use super::{iter_next, Element, ElementId};
use crate::{element_id, mesh::Mesh, INVALID_IND};

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
        Self { id, mesh, phantom: PhantomData}
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
