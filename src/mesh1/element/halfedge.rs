use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    mesh1::{element::{EdgeId, VertexId}, mesh::MeshCore}, INVALID_IND
};

use super::{ElementId, FaceId};

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
    fn edge(&self) -> EdgeId;
    fn sibling(&self) -> HalfedgeId;
    fn incoming_next(&self) -> HalfedgeId;
}

pub struct HalfedgeIter<'m, M: MeshCore> {
    pub id: HalfedgeId,
    pub data: &'m M::Halfedge,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m M>,
}

pub struct HalfedgeIterMut<'m, M: MeshCore> {
    pub id: HalfedgeId,
    pub data: &'m M::Halfedge,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m mut M>,
}

impl<'m, M: MeshCore> HalfedgeIter<'m, M> {
    pub fn new(id: HalfedgeId, mesh: &'m M) -> Self {
        let data = mesh.halfedge(id);
        HalfedgeIter {
            id,
            data,
            mesh: NonNull::from(mesh),
            _marker: PhantomData,
        }
    }
}

impl<'m, M: MeshCore> HalfedgeIterMut<'m, M> {
    pub fn new(id: HalfedgeId, mesh: &'m mut M) -> Self {
        let mut mesh = NonNull::from_mut(mesh);
        unsafe {
            let data = mesh.as_mut().halfedge_mut(id);
            HalfedgeIterMut {
                id,
                data,
                mesh,
                _marker: PhantomData,
            }
        }
    }
}
