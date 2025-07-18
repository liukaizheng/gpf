use std::{iter::from_fn, ops::Deref};

use crate::{
    INVALID_IND,
    mesh1::{
        element::HalfedgeIter,
        mesh::{Mesh, MeshCore},
    },
};

use super::{ElementId, Halfedge, HalfedgeExt, HalfedgeId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(usize);

impl Default for VertexId {
    #[inline]
    fn default() -> Self {
        VertexId(INVALID_IND)
    }
}

impl From<usize> for VertexId {
    #[inline]
    fn from(index: usize) -> Self {
        VertexId(index)
    }
}

impl ElementId for VertexId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for VertexId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Vertex {
    fn halfedge(&self) -> HalfedgeId;
    fn set_halfedge(&mut self, hid: HalfedgeId);
}

pub struct VertexIter<'m, M: MeshCore> {
    id: VertexId,
    data: &'m M::Vertex,
    mesh: &'m M,
}

impl<'m, M: MeshCore> VertexIter<'m, M> {
    pub fn new(id: VertexId, mesh: &'m M) -> Self {
        let data = mesh.vertex(id);
        VertexIter { id, data, mesh }
    }
}

pub struct CirculateVertex<'m, M: Mesh> {
    id: VertexId,
    first_hid: HalfedgeId,
    hid: HalfedgeId,
    he: &'m M::Halfedge,
    mesh: &'m M,
    first: bool,
}

trait HalfedgeIncomingNext {
    fn next(&mut self);
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> CirculateVertex<'m, M> {
    #[inline]
    pub fn new(id: VertexId, from_hid: HalfedgeId, mesh: &'m M) -> Self {
        let hid = mesh.halfedge(from_hid).prev();
        CirculateVertex {
            id,
            first: true,
            first_hid: hid,
            hid,
            he: &mesh.halfedge(hid),
            mesh,
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for CirculateVertex<'m, M> {
    default fn next(&mut self) {
        self.hid = self.mesh.he_incoming_next(self.hid);
        self.he = self.mesh.halfedge(self.hid);
        self.first = false;
    }
}

impl<'m, M: Mesh<Halfedge: HalfedgeExt>> HalfedgeIncomingNext for CirculateVertex<'m, M> {
    fn next(&mut self) {
        self.hid = self.he.incoming_next();
        self.he = self.mesh.halfedge(self.hid);
        self.first = false;
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexIter<'m, M> {
    pub fn incoming_halfedges(&self) -> impl Iterator<Item = HalfedgeIter<'m, M>> {
        let mut cv = CirculateVertex::new(self.id, self.data.halfedge(), self.mesh);
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = HalfedgeIter {
                id: cv.hid,
                mesh: self.mesh,
            };
            cv.next();
            Some(ret)
        })
    }
    pub fn incoming_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = CirculateVertex::new(self.id, self.data.halfedge(), self.mesh);
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = cv.hid;
            debug_assert!(cv.he.vertex() == self.id);
            cv.next();
            Some(ret)
        })
    }

    pub fn outgoing_halfedges(&self) -> impl Iterator<Item = HalfedgeIter<'m, M>> {
        let mut cv = CirculateVertex::new(self.id, self.data.halfedge(), self.mesh);
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = HalfedgeIter {
                id: cv.he.next(),
                mesh: self.mesh,
            };
            cv.next();
            Some(ret)
        })
    }

    pub fn outgoing_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = CirculateVertex::new(self.id, self.data.halfedge(), self.mesh);
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = cv.he.next();
            debug_assert!(cv.he.vertex() == self.id);
            cv.next();
            Some(ret)
        })
    }
}
