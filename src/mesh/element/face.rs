use std::ops::{Add, Deref, DerefMut, Index, IndexMut, Mul};

use super::{
    Element, ElementId, ElementIndex, Halfedge, HalfedgeId, Vertex, WHIter, Wire, iter_next,
};
use crate::{
    INVALID_IND, element_id,
    mesh::{HoleAwareMesh, Mesh},
};

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

pub struct FHIterImpl<'a, M: Mesh> {
    first_loop_hid: HalfedgeId,
    current: WHIter<'a, M>,
    first: bool,
}

impl<'a, M: Mesh> FHIterImpl<'a, M> {
    pub fn new(mesh: &'a M, fid: FaceId) -> Self {
        let first_loop_hid = mesh.f_halfedge(fid);
        Self {
            first_loop_hid,
            current: WHIter::new(mesh, first_loop_hid),
            first: true,
        }
    }

    #[inline]
    fn prev(&mut self) {
        self.current.0.prev();
        if self.current.is_end() {
            let mesh = self.current.0.mesh;
            let first_hid = self.current.0.first_hid;
            self.current = WHIter::new(mesh, mesh.f_loop_next_first_halfedge(first_hid));
            self.first = false;
        }
    }
}

impl<'a, M: Mesh> Element for FHIterImpl<'a, M> {
    type Item = HalfedgeId;

    #[inline]
    fn item(&self) -> Self::Item {
        self.current.0.current
    }

    #[inline]
    fn valid(&self) -> bool {
        true
    }

    #[inline]
    fn next(&mut self) {
        Element::next(&mut self.current);
        if self.current.is_end() {
            let mesh = self.current.0.mesh;
            let first_hid = self.current.0.first_hid;
            self.current = WHIter::new(mesh, mesh.f_loop_next_first_halfedge(first_hid));
            self.first = false;
        }
    }

    fn is_end(&self) -> bool {
        !self.first && self.current.0.first && self.current.0.current == self.first_loop_hid
    }
}

pub struct FHIter<'a, M: Mesh>(FHIterImpl<'a, M>);

impl<'a, M: Mesh> Element for FHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.0.current.0.mesh, self.0.current.0.current)
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
            self.0.current.0.mesh,
            self.0.current.0.mesh.he_to(self.0.current.0.current),
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

impl<'a, A: Allocator + Copy> Face<'a, HoleAwareMesh<A>> {
    pub fn wires(&self) -> FaceWireIter<'a, A> {
        FaceWireIter::new(self.mesh, self.id)
    }
}

pub struct FaceWireIter<'a, A: Allocator + Copy> {
    mesh: &'a HoleAwareMesh<A>,
    first_hid: HalfedgeId,
    curr_hid: HalfedgeId,
    first: bool,
}

impl<'a, A: Allocator + Copy> FaceWireIter<'a, A> {
    fn new(mesh: &'a HoleAwareMesh<A>, face_id: FaceId) -> Self {
        let first_hid = mesh.f_halfedge(face_id);
        Self {
            mesh,
            first_hid,
            curr_hid: first_hid,
            first: true,
        }
    }
}

impl<'a, A: Allocator + Copy> Element for FaceWireIter<'a, A> {
    type Item = Wire<'a, HoleAwareMesh<A>>;

    #[inline(always)]
    fn item(&self) -> Self::Item {
        Wire::new(self.mesh, self.curr_hid)
    }

    #[inline(always)]
    fn valid(&self) -> bool {
        true
    }

    #[inline(always)]
    fn next(&mut self) {
        self.curr_hid = self.mesh.f_loop_next_first_halfedge(self.curr_hid);
    }

    #[inline(always)]
    fn is_end(&self) -> bool {
        self.curr_hid == self.first_hid && !self.first
    }
}

impl<'a, A: Allocator + Copy> Iterator for FaceWireIter<'a, A> {
    type Item = Wire<'a, HoleAwareMesh<A>>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
