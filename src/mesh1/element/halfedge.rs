use std::ops::Deref;

use crate::{INVALID_IND, mesh1::mesh::Mesh};

use super::{ElementId, FaceId, VertexId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(pub usize);

impl Default for HalfedgeId {
    #[inline]
    fn default() -> Self {
        HalfedgeId(INVALID_IND)
    }
}

impl From<usize> for HalfedgeId {
    #[inline]
    fn from(index: usize) -> Self {
        HalfedgeId(index)
    }
}

impl ElementId for HalfedgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for HalfedgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Halfedge {
    fn vertex(&self) -> VertexId;
    fn set_vertex(&mut self, vertex: VertexId);

    fn prev(&self) -> HalfedgeId;
    fn set_prev(&mut self, prev: HalfedgeId);

    fn next(&self) -> HalfedgeId;
    fn set_next(&mut self, next: HalfedgeId);

    fn face(&self) -> FaceId;
    fn set_face(&mut self, face: FaceId);
}

pub trait HalfedgeExt: Halfedge {
    fn sibling(&self) -> HalfedgeId;
    fn incoming_next(&self) -> HalfedgeId;
}

pub struct HalfedgeIter<'m, M: Mesh> {
    pub id: HalfedgeId,
    pub mesh: &'m M,
}
