use std::{iter::from_fn, marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND,
    mesh1::{
        element::{EdgeId, EdgeIter, EdgeIterMut, HalfedgeIter, HalfedgeIterMut},
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
    pub id: VertexId,
    pub data: &'m M::Vertex,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m M>,
}

pub struct VertexIterMut<'m, M: MeshCore> {
    pub id: VertexId,
    pub data: &'m M::Vertex,
    pub mesh: NonNull<M>,
    _marker: PhantomData<&'m mut M>,
}

impl<'m, M: MeshCore> VertexIter<'m, M> {
    pub fn new(id: VertexId, mesh: &'m M) -> Self {
        let data = mesh.vertex(id);
        VertexIter {
            id,
            data,
            mesh: NonNull::from(mesh),
            _marker: PhantomData,
        }
    }
}

impl<'m, M: MeshCore> VertexIterMut<'m, M> {
    pub fn new(id: VertexId, mesh: &'m mut M) -> Self {
        let mut mesh = NonNull::from_mut(mesh);
        unsafe {
            let data = mesh.as_mut().vertex_mut(id);
            VertexIterMut {
                id,
                data,
                mesh,
                _marker: PhantomData,
            }
        }
    }
}

pub struct CirculateVertex<'m, M: Mesh> {
    first_hid: HalfedgeId,
    hid: HalfedgeId,
    he: &'m M::Halfedge,
    mesh: NonNull<M>,
    first: bool,
    _marker: PhantomData<&'m M>,
}

pub struct CirculateVertexMut<'m, M: Mesh> {
    first_hid: HalfedgeId,
    hid: HalfedgeId,
    he: &'m mut M::Halfedge,
    mesh: NonNull<M>,
    first: bool,
    _marker: PhantomData<&'m mut M>,
}

trait HalfedgeIncomingNext {
    fn next(&mut self);
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> CirculateVertex<'m, M> {
    #[inline]
    pub fn new(from_hid: HalfedgeId, mesh: &M) -> Self {
        let mesh = NonNull::from(mesh);
        unsafe {
            let hid = mesh.as_ref().halfedge(from_hid).prev();
            CirculateVertex {
                first: true,
                first_hid: hid,
                hid,
                he: &mesh.as_ref().halfedge(hid),
                mesh,
                _marker: PhantomData,
            }
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for CirculateVertex<'m, M> {
    default fn next(&mut self) {
        let mesh = unsafe { self.mesh.as_ref() };
        self.hid = mesh.he_incoming_next(self.hid);
        self.he = mesh.halfedge(self.hid);
        self.first = false;
    }
}

impl<'m, M: Mesh<Halfedge: HalfedgeExt>> HalfedgeIncomingNext for CirculateVertex<'m, M> {
    fn next(&mut self) {
        self.hid = self.he.incoming_next();
        unsafe {
            self.he = self.mesh.as_ref().halfedge(self.hid);
        }
        self.first = false;
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> CirculateVertexMut<'m, M> {
    #[inline]
    pub fn new(from_hid: HalfedgeId, mesh: &'m mut M) -> Self {
        let mut mesh = NonNull::from_mut(mesh);
        unsafe {
            let hid = mesh.as_ref().halfedge(from_hid).prev();
            CirculateVertexMut {
                first: true,
                first_hid: hid,
                hid,
                he: mesh.as_mut().halfedge_mut(hid),
                mesh,
                _marker: PhantomData,
            }
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for CirculateVertexMut<'m, M> {
    default fn next(&mut self) {
        unsafe {
            self.hid = self.mesh.as_ref().he_incoming_next(self.hid);
            self.he = self.mesh.as_mut().halfedge_mut(self.hid);
            self.first = false;
        }
    }
}

impl<'m, M: Mesh<Halfedge: HalfedgeExt>> HalfedgeIncomingNext for CirculateVertexMut<'m, M> {
    fn next(&mut self) {
        self.hid = self.he.incoming_next();
        unsafe {
            self.he = self.mesh.as_mut().halfedge_mut(self.hid);
        }
        self.first = false;
    }
}

pub trait VertexEdgesAndVertices<'m, M: Mesh<Edge: 'm> + 'm> {
    fn edges(&self) -> impl Iterator<Item = EdgeIter<'m, M>>;
    fn edge_ids(&self) -> impl Iterator<Item = EdgeId>;
    fn vertices(&self) -> impl Iterator<Item = VertexIter<'m, M>>;
    fn vertex_ids(&self) -> impl Iterator<Item = VertexId>;
}

pub trait VertexEdgesAndVerticesMut<'m, M: Mesh + 'm>: VertexEdgesAndVertices<'m, M> {
    fn edges_mut(&mut self) -> impl Iterator<Item = EdgeIterMut<'m, M>>;
    fn vertices_mut(&mut self) -> impl Iterator<Item = VertexIterMut<'m, M>>;
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexIter<'m, M> {
    pub fn incoming_halfedges(&self) -> impl Iterator<Item = HalfedgeIter<'m, M>> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIter::new(cv.hid, self.mesh.as_ref()) };
            cv.next();
            Some(ret)
        })
    }
    pub fn incoming_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
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
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIter::new(cv.he.next(), self.mesh.as_ref()) };
            cv.next();
            Some(ret)
        })
    }

    pub fn outgoing_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
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

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexIterMut<'m, M> {
    pub fn incoming_halfedges(&self) -> impl Iterator<Item = HalfedgeIter<'m, M>> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIter::new(cv.hid, self.mesh.as_ref()) };
            cv.next();
            Some(ret)
        })
    }
    pub fn incoming_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
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
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIter::new(cv.he.next(), self.mesh.as_ref()) };
            cv.next();
            Some(ret)
        })
    }

    pub fn outgoing_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
        let mut cv = unsafe { CirculateVertex::new(self.data.halfedge(), self.mesh.as_ref()) };
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
    pub fn halfedge_mut(&mut self) -> HalfedgeIterMut<'m, M> {
        unsafe { HalfedgeIterMut::new(self.data.halfedge(), self.mesh.as_mut()) }
    }

    pub fn incoming_halfedges_mut(&mut self) -> impl Iterator<Item = HalfedgeIterMut<'m, M>> {
        let mut cv = unsafe { CirculateVertexMut::new(self.data.halfedge(), self.mesh.as_mut()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIterMut::new(cv.hid, self.mesh.as_mut()) };
            cv.next();
            Some(ret)
        })
    }

    pub fn outgoing_halfedges_mut(&mut self) -> impl Iterator<Item = HalfedgeIterMut<'m, M>> {
        let mut cv = unsafe { CirculateVertexMut::new(self.data.halfedge(), self.mesh.as_mut()) };
        from_fn(move || {
            if !cv.valid() {
                return None;
            }
            let ret = unsafe { HalfedgeIterMut::new(cv.he.next(), self.mesh.as_mut()) };
            cv.next();
            Some(ret)
        })
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexEdgesAndVertices<'m, M>
    for VertexIter<'m, M>
{
    fn edges(&self) -> impl Iterator<Item = EdgeIter<'m, M>> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(EdgeIter::new(eid, mesh));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                prev = Some(he_item);
                                return Some(EdgeIter::new(eid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(eid);
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                prev = Some(he_item);
                                return Some(eid);
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn vertices(&self) -> impl Iterator<Item = VertexIter<'m, M>> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(VertexIter::new(mesh.he_from(hid), mesh));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                let vid = he_item.data.vertex();
                                prev = Some(he_item);
                                return Some(VertexIter::new(vid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn vertex_ids(&self) -> impl Iterator<Item = VertexId> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(mesh.he_from(hid));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                let vid = he_item.data.vertex();
                                prev = Some(he_item);
                                return Some(vid);
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexEdgesAndVertices<'m, M>
    for VertexIterMut<'m, M>
{
    fn edges(&self) -> impl Iterator<Item = EdgeIter<'m, M>> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(EdgeIter::new(eid, mesh));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                prev = Some(he_item);
                                return Some(EdgeIter::new(eid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(eid);
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                prev = Some(he_item);
                                return Some(eid);
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn vertices(&self) -> impl Iterator<Item = VertexIter<'m, M>> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(VertexIter::new(mesh.he_from(hid), mesh));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                let vid = he_item.data.vertex();
                                prev = Some(he_item);
                                return Some(VertexIter::new(vid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn vertex_ids(&self) -> impl Iterator<Item = VertexId> {
        let mut he_iter = self.incoming_halfedges();
        let mut prev: Option<HalfedgeIter<'m, M>> = None;
        unsafe {
            let mesh = self.mesh.as_ref();
            from_fn(move || {
                loop {
                    if let Some(he_item) = prev.take() {
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(mesh.he_from(hid));
                        }
                    } else {
                        if let Some(he_item) = he_iter.next() {
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                let vid = he_item.data.vertex();
                                prev = Some(he_item);
                                return Some(vid);
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexEdgesAndVerticesMut<'m, M>
    for VertexIterMut<'m, M>
{
    fn edges_mut(&mut self) -> impl Iterator<Item = EdgeIterMut<'m, M>> {
        let mut he_iter = self.incoming_halfedges_mut();
        let mut prev: Option<HalfedgeIterMut<'m, M>> = None;
        unsafe {
            from_fn(move || {
                loop {
                    if let Some(mut he_item) = prev.take() {
                        let mesh = he_item.mesh.as_mut();
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(EdgeIterMut::new(eid, mesh));
                        }
                    } else {
                        if let Some(mut he_item) = he_iter.next() {
                            let mesh = he_item.mesh.as_mut();
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                prev = Some(he_item);
                                return Some(EdgeIterMut::new(eid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }

    fn vertices_mut(&mut self) -> impl Iterator<Item = VertexIterMut<'m, M>> {
        let mut he_iter = self.incoming_halfedges_mut();
        let mut prev: Option<HalfedgeIterMut<'m, M>> = None;
        unsafe {
            from_fn(move || {
                loop {
                    if let Some(mut he_item) = prev.take() {
                        let mesh = he_item.mesh.as_mut();
                        let hid = he_item.data.next();
                        let eid = mesh.he_edge(hid);
                        if mesh.e_halfedge(eid) == hid {
                            return Some(VertexIterMut::new(mesh.he_from(hid), mesh));
                        }
                    } else {
                        if let Some(mut he_item) = he_iter.next() {
                            let mesh = he_item.mesh.as_mut();
                            let eid = mesh.he_edge(he_item.id);
                            if mesh.e_halfedge(eid) == he_item.id {
                                let vid = he_item.data.vertex();
                                prev = Some(he_item);
                                return Some(VertexIterMut::new(vid, mesh));
                            } else {
                                prev = Some(he_item);
                            }
                        } else {
                            return None;
                        }
                    }
                }
            })
        }
    }
}
