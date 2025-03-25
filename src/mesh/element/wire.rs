use crate::mesh::Mesh;

use super::{Element, Halfedge, HalfedgeId, Vertex, iter_next};

pub struct Wire<'a, M: Mesh> {
    mesh: &'a M,
    hid: HalfedgeId,
}

impl<'a, M: Mesh> Wire<'a, M> {
    pub fn new(mesh: &'a M, hid: HalfedgeId) -> Self {
        Wire { mesh, hid }
    }

    pub fn halfedges(&self) -> WHIter<'a, M> {
        WHIter::new(self.mesh, self.hid)
    }

    pub fn vertices(&self) -> WVIter<'a, M> {
        WVIter::new(self.mesh, self.hid)
    }
}

/// Iterator over the halfedges of a face wire.
pub struct WHIterImpl<'a, M: Mesh> {
    pub(crate) mesh: &'a M,
    pub(crate) first_hid: HalfedgeId,
    pub(crate) current: HalfedgeId,
    pub(crate) first: bool,
}

impl<'a, M: Mesh> WHIterImpl<'a, M> {
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

impl<'a, M: Mesh> Element for WHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline]
    fn item(&self) -> Self::Item {
        Halfedge::new(self.0.mesh, self.0.current)
    }

    #[inline]
    fn valid(&self) -> bool {
        true
    }

    #[inline]
    fn next(&mut self) {
        self.0.next();
    }

    #[inline]
    fn is_end(&self) -> bool {
        self.0.is_end()
    }
}
pub struct WHIter<'a, M: Mesh>(pub(crate) WHIterImpl<'a, M>);

impl<'a, M: Mesh> WHIter<'a, M> {
    pub fn new(mesh: &'a M, hid: HalfedgeId) -> Self {
        Self(WHIterImpl::new(mesh, hid))
    }
}

impl<'a, M: Mesh> Iterator for WHIter<'a, M> {
    type Item = Halfedge<'a, M>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
impl<'a, M: Mesh> DoubleEndedIterator for WHIter<'a, M> {
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

pub struct WVIter<'a, M: Mesh>(WHIterImpl<'a, M>);
impl<'a, M: Mesh> WVIter<'a, M> {
    pub fn new(mesh: &'a M, hid: HalfedgeId) -> Self {
        Self(WHIterImpl::new(mesh, hid))
    }
}

impl<'a, M: Mesh> Element for WVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline]
    fn item(&self) -> Self::Item {
        Vertex::new(self.0.mesh, self.0.mesh.he_to(self.0.current))
    }

    #[inline]
    fn valid(&self) -> bool {
        true
    }

    #[inline]
    fn next(&mut self) {
        self.0.next();
    }

    #[inline]
    fn is_end(&self) -> bool {
        self.0.is_end()
    }
}

impl<'a, M: Mesh> Iterator for WVIter<'a, M> {
    type Item = Vertex<'a, M>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
impl<'a, M: Mesh> DoubleEndedIterator for WVIter<'a, M> {
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
