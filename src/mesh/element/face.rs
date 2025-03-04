use std::ops::{Add, Deref, DerefMut, Index, IndexMut, Mul};

use super::{iter_next, Element, ElementId, ElementIndex, Halfedge, HalfedgeId, Vertex};
use crate::{element_id, mesh::Mesh, INVALID_IND};

use std::alloc::Allocator;

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
        FHIter(FHIterImpl::new(self.mesh, self.id))
    }

    pub fn vertices(&self) -> FVIter<'a, M> {
        FVIter(FHIterImpl::new(self.mesh, self.id))
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
        self.mesh.f_is_valid(self.id)
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

/// Iterator over the halfedges of a face loop.
pub struct LHIter<'a, M: Mesh> {
    mesh: &'a M,
    first_hid: HalfedgeId,
    current: HalfedgeId,
    first: bool,
}

impl<'a, M: Mesh> LHIter<'a, M> {
    pub fn new(mesh: &'a M, hid: HalfedgeId) -> Self {
        Self {
            mesh,
            first_hid: hid,
            current: hid,
            first: true,
        }
    }

    #[inline]
    pub fn prev(&mut self) {
        self.current = self.mesh.he_prev(self.current);
        self.first = false;
    }
}

impl<'a, M: Mesh> Element for LHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.mesh, self.current)
    }

    #[inline]
    fn valid(&self) -> bool {
        true
    }

    #[inline]
    fn next(&mut self) {
        self.current = self.mesh.he_next(self.current);
        self.first = false;
    }

    #[inline]
    fn is_end(&self) -> bool {
        self.current == self.first_hid && !self.first
    }
}

impl<'a, M: Mesh> Iterator for LHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
impl<'a, M: Mesh> DoubleEndedIterator for LHIter<'a, M> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.is_end() {
            return None;
        } else {
            let ret = Some(self.item());
            self.prev();
            ret
        }
    }
}

pub struct FHIterImpl<'a, M: Mesh> {
    first_loop_hid: HalfedgeId,
    current: LHIter<'a, M>,
    first: bool,
}

impl<'a, M: Mesh> FHIterImpl<'a, M> {
    pub fn new(mesh: &'a M, fid: FaceId) -> Self {
        let first_loop_hid = mesh.f_halfedge(fid);
        Self {
            first_loop_hid,
            current: LHIter::new(mesh, first_loop_hid),
            first: true,
        }
    }

    #[inline]
    fn prev(&mut self) {
        self.current.prev();
        if self.current.is_end() {
            let mesh = self.current.mesh;
            let first_hid = self.current.first_hid;
            self.current = LHIter::new(mesh, mesh.f_loop_next_first_halfedge(first_hid));
            self.first = false;
        }
    }
}

impl<'a, M: Mesh> Element for FHIterImpl<'a, M> {
    type Item = HalfedgeId;

    #[inline]
    fn item(&self) -> Self::Item {
        self.current.current
    }

    #[inline]
    fn valid(&self) -> bool {
        true
    }

    #[inline]
    fn next(&mut self) {
        Element::next(&mut self.current);
        if self.current.is_end() {
            let mesh = self.current.mesh;
            let first_hid = self.current.first_hid;
            self.current = LHIter::new(mesh, mesh.f_loop_next_first_halfedge(first_hid));
            self.first = false;
        }
    }

    fn is_end(&self) -> bool {
        !self.first && self.current.first && self.current.current == self.first_loop_hid
    }
}

pub struct FHIter<'a, M: Mesh>(FHIterImpl<'a, M>);

impl<'a, M: Mesh> Element for FHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.0.current.mesh, self.0.current.current)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.0.next();
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.0.is_end()
    }
}

impl<'a, M: Mesh> Iterator for FHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

impl<'a, M: Mesh> DoubleEndedIterator for FHIter<'a, M> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.is_end() {
            return None;
        } else {
            let ret = Some(self.item());
            self.0.prev();
            ret
        }
    }
}

pub struct FVIter<'a, M: Mesh>(FHIterImpl<'a, M>);

impl<'a, M: Mesh> Element for FVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Vertex::new(
            self.0.current.mesh,
            self.0.current.mesh.he_to(self.0.current.current),
        )
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.0.next();
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.0.is_end()
    }
}

impl<'a, M: Mesh> Iterator for FVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}

impl<'a, M: Mesh> DoubleEndedIterator for FVIter<'a, M> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.is_end() {
            return None;
        } else {
            let ret = Some(self.item());
            self.0.prev();
            ret
        }
    }
}
