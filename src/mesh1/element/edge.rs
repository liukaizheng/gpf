use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{HalfedgeId, HalfedgeIter},
        mesh::Mesh,
    },
};

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

element_iter_struct!(struct EdgeIter -> Mesh, EdgeId, Edge, from, as_ref, edge, {});
element_iter_struct!(struct EdgeIterMut -> Mesh, EdgeId, Edge, from_mut, as_mut, edge_mut, {mut});

struct EdgeHalfedges<'m, M: Mesh> {
    first_hid: HalfedgeId,
    hid: HalfedgeId,
    he: &'m M::Halfedge,
    mesh: NonNull<M>,
    first: bool,
    _marker: PhantomData<&'m M>,
}

impl<'m, M: Mesh> EdgeHalfedges<'m, M> {
    #[inline]
    fn new(first_hid: HalfedgeId, mesh: &'m M) -> Self {
        Self {
            first_hid,
            hid: first_hid,
            he: mesh.halfedge(first_hid),
            mesh: NonNull::from(mesh),
            first: true,
            _marker: PhantomData,
        }
    }

    #[inline]
    fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh> Iterator for EdgeHalfedges<'m, M> {
    type Item = HalfedgeIter<'m, M>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if !self.valid() {
            return None;
        }
        self.first = false;
        unsafe {
            let ret = HalfedgeIter::new(self.hid, self.mesh.as_ref());
            self.hid = self.mesh.as_ref().he_next(self.hid);
            Some(ret)
        }
    }
}

pub trait EdgeMethod<'m, M: Mesh> {
    fn halfedge(&self) -> HalfedgeIter<'m, M>;
    fn halfedges(&self) -> IntoIter<HalfedgeIter<'m, M>>;
}
