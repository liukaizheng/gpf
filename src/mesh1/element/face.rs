use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND,
    mesh1::mesh::{Mesh, MeshCore},
};

use super::{ElementId, Halfedge, HalfedgeId};

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

pub struct FaceIter<'m, M: MeshCore> {
    pub id: FaceId,
    pub data: &'m M::Face,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m M>,
}

pub struct FaceIterMut<'m, M: MeshCore> {
    pub id: FaceId,
    pub data: &'m M::Face,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m mut M>,
}

impl<'m, M: MeshCore> FaceIter<'m, M> {
    pub fn new(id: FaceId, mesh: &'m M) -> Self {
        let data = mesh.face(id);
        FaceIter {
            id,
            data,
            mesh: NonNull::from(mesh),
            _marker: PhantomData,
        }
    }
}

impl<'m, M: MeshCore> FaceIterMut<'m, M> {
    pub fn new(id: FaceId, mesh: &'m mut M) -> Self {
        let mut mesh = NonNull::from_mut(mesh);
        unsafe {
            let data = mesh.as_mut().face_mut(id);
            FaceIterMut {
                id,
                data,
                mesh,
                _marker: PhantomData,
            }
        }
    }
}
