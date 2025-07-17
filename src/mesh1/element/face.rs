use std::ops::Deref;

use crate::{mesh1::mesh::Mesh, INVALID_IND};

use super::{ ElementId, Halfedge, HalfedgeId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

impl Default for FaceId {
    #[inline]
    fn default() -> Self {
        FaceId(INVALID_IND)
    }
}

impl From<usize> for FaceId {
    #[inline]
    fn from(index: usize) -> Self {
        FaceId(index)
    }
}

impl ElementId for FaceId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for FaceId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Face {
    fn halfedge(&self) -> HalfedgeId;
    fn set_halfedge(&mut self, halfedge: HalfedgeId);
}

pub struct FaceIter<'m, M: Mesh> {
    id: FaceId,
    mesh:&'m M,
}
