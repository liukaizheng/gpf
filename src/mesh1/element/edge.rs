
use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{mesh1::{element::HalfedgeId, mesh::Mesh}, INVALID_IND};

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

impl Default for EdgeId {
    #[inline]
    fn default() -> Self {
        EdgeId(INVALID_IND)
    }
}

impl From<usize> for EdgeId {
    #[inline]
    fn from(index: usize) -> Self {
        EdgeId(index)
    }
}

impl ElementId for EdgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for EdgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Edge {
    fn halfedge(&self) -> HalfedgeId;
}

pub struct EdgeIter<'m, M: Mesh> {
    pub id: EdgeId,
    pub data: &'m M::Edge,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m M>,
}

pub struct EdgeIterMut<'m, M: Mesh> {
    pub id: EdgeId,
    pub data: &'m M::Edge,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m mut M>,
}

impl<'m, M: Mesh> EdgeIter<'m, M> {
    pub fn new(id: EdgeId, mesh: &'m M) -> Self {
        let data = mesh.edge(id);
        EdgeIter {
            id,
            data,
            mesh: NonNull::from(mesh),
            _marker: PhantomData,
        }
    }
}

impl<'m, M: Mesh> EdgeIterMut<'m, M> {
    pub fn new(id: EdgeId, mesh: &'m mut M) -> Self {
        let mut mesh = NonNull::from_mut(mesh);
        unsafe {
            let data = mesh.as_mut().edge_mut(id);
            EdgeIterMut {
                id,
                data,
                mesh,
                _marker: PhantomData,
            }
        }
    }
}
