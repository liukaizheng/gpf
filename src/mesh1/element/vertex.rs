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

// Macro to generate vertex iterator structs
macro_rules! vertex_iter_struct {
    ($name:ident, $marker:ty) => {
        pub struct $name<'m, M: MeshCore> {
            pub id: VertexId,
            pub data: &'m M::Vertex,
            pub mesh: NonNull<M>,
            _marker: PhantomData<$marker>,
        }
    };
}

vertex_iter_struct!(VertexIter, &'m M);
vertex_iter_struct!(VertexIterMut, &'m mut M);

// Macro to generate circulate vertex structs
macro_rules! circulate_vertex_struct {
    ($name:ident, $he_type:ty, $marker:ty) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            hid: HalfedgeId,
            he: $he_type,
            mesh: NonNull<M>,
            first: bool,
            _marker: PhantomData<$marker>,
        }
    };
}

circulate_vertex_struct!(CirculateVertex, &'m M::Halfedge, &'m M);
circulate_vertex_struct!(CirculateVertexMut, &'m mut M::Halfedge, &'m mut M);

trait HalfedgeIncomingNext {
    fn next(&mut self);
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
        let mesh_ptr = NonNull::from(&*mesh);
        let data = mesh.vertex(id);
        VertexIterMut {
            id,
            data,
            mesh: mesh_ptr,
            _marker: PhantomData,
        }
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> CirculateVertex<'m, M> {
    #[inline]
    pub fn new(from_hid: HalfedgeId, mesh: &'m M) -> Self {
        let mesh_ptr = NonNull::from(mesh);
        let hid = mesh.halfedge(from_hid).prev();
        let he = mesh.halfedge(hid);
        CirculateVertex {
            first: true,
            first_hid: hid,
            hid,
            he,
            mesh: mesh_ptr,
            _marker: PhantomData,
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> CirculateVertexMut<'m, M> {
    #[inline]
    pub fn new(from_hid: HalfedgeId, mesh: &'m mut M) -> Self {
        let mesh_ptr = NonNull::from(&mut *mesh);
        let hid = mesh.halfedge(from_hid).prev();
        let he = mesh.halfedge_mut(hid);
        CirculateVertexMut {
            first: true,
            first_hid: hid,
            hid,
            he,
            mesh: mesh_ptr,
            _marker: PhantomData,
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

// Macro to generate HalfedgeIncomingNext implementations
macro_rules! halfedge_incoming_next_impl {
    ($name:ident, $he_access:ident) => {
        impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for $name<'m, M> {
            default fn next(&mut self) {
                unsafe {
                    let mesh = self.mesh.as_ref();
                    self.hid = mesh.he_incoming_next(self.hid);
                    self.he = mesh.$he_access(self.hid);
                    self.first = false;
                }
            }
        }

        impl<'m, M: Mesh<Halfedge: HalfedgeExt>> HalfedgeIncomingNext for $name<'m, M> {
            fn next(&mut self) {
                self.hid = self.he.incoming_next();
                unsafe {
                    self.he = self.mesh.as_ref().$he_access(self.hid);
                }
                self.first = false;
            }
        }
    };
}

halfedge_incoming_next_impl!(CirculateVertex, halfedge);

// Special implementation for CirculateVertexMut since it needs mutable access
impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for CirculateVertexMut<'m, M> {
    default fn next(&mut self) {
        unsafe {
            let mesh = self.mesh.as_ref();
            self.hid = mesh.he_incoming_next(self.hid);
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

// Macro to generate halfedge iteration methods
macro_rules! halfedge_iter_methods {
    ($impl_type:ident, $circulate_type:ident, $he_iter_type:ident) => {
        impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> $impl_type<'m, M> {
            pub fn incoming_halfedges(&self) -> impl Iterator<Item = $he_iter_type<'m, M>> {
                let mut cv =
                    unsafe { $circulate_type::new(self.data.halfedge(), self.mesh.as_ref()) };
                from_fn(move || {
                    if !cv.valid() {
                        return None;
                    }
                    let ret = unsafe { $he_iter_type::new(cv.hid, self.mesh.as_ref()) };
                    cv.next();
                    Some(ret)
                })
            }

            pub fn incoming_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
                let mut cv =
                    unsafe { $circulate_type::new(self.data.halfedge(), self.mesh.as_ref()) };
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

            pub fn outgoing_halfedges(&self) -> impl Iterator<Item = $he_iter_type<'m, M>> {
                let mut cv =
                    unsafe { $circulate_type::new(self.data.halfedge(), self.mesh.as_ref()) };
                from_fn(move || {
                    if !cv.valid() {
                        return None;
                    }
                    let ret = unsafe { $he_iter_type::new(cv.he.next(), self.mesh.as_ref()) };
                    cv.next();
                    Some(ret)
                })
            }

            pub fn outgoing_halfedge_ids(&self) -> impl Iterator<Item = HalfedgeId> {
                let mut cv =
                    unsafe { $circulate_type::new(self.data.halfedge(), self.mesh.as_ref()) };
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
    };
}

halfedge_iter_methods!(VertexIter, CirculateVertex, HalfedgeIter);
halfedge_iter_methods!(VertexIterMut, CirculateVertex, HalfedgeIter);

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexIterMut<'m, M> {
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

// Macro to generate edge/vertex iteration logic
macro_rules! edge_vertex_iter_impl {
    ($for_type:ident, $edge_iter_type:ident, $vertex_iter_type:ident) => {
        impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexEdgesAndVertices<'m, M>
            for $for_type<'m, M>
        {
            fn edges(&self) -> impl Iterator<Item = $edge_iter_type<'m, M>> {
                let mut he_iter = self.incoming_halfedges();
                let mut prev: Option<HalfedgeIter<'m, M>> = None;
                let mesh = unsafe { self.mesh.as_ref() };
                from_fn(move || {
                    loop {
                        if let Some(he_item) = prev.take() {
                            let hid = he_item.data.next();
                            let eid = mesh.he_edge(hid);
                            if mesh.e_halfedge(eid) == hid {
                                return Some($edge_iter_type::new(eid, mesh));
                            }
                        } else {
                            if let Some(he_item) = he_iter.next() {
                                let eid = mesh.he_edge(he_item.id);
                                if mesh.e_halfedge(eid) == he_item.id {
                                    prev = Some(he_item);
                                    return Some($edge_iter_type::new(eid, mesh));
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

            fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
                let mut he_iter = self.incoming_halfedges();
                let mut prev: Option<HalfedgeIter<'m, M>> = None;
                let mesh = unsafe { self.mesh.as_ref() };
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

            fn vertices(&self) -> impl Iterator<Item = $vertex_iter_type<'m, M>> {
                let mut he_iter = self.incoming_halfedges();
                let mut prev: Option<HalfedgeIter<'m, M>> = None;
                let mesh = unsafe { self.mesh.as_ref() };
                from_fn(move || {
                    loop {
                        if let Some(he_item) = prev.take() {
                            let hid = he_item.data.next();
                            let eid = mesh.he_edge(hid);
                            if mesh.e_halfedge(eid) == hid {
                                return Some($vertex_iter_type::new(mesh.he_to(hid), mesh));
                            }
                        } else {
                            if let Some(he_item) = he_iter.next() {
                                let eid = mesh.he_edge(he_item.id);
                                if mesh.e_halfedge(eid) == he_item.id {
                                    let vid = mesh.he_to(he_item.data.prev());
                                    prev = Some(he_item);
                                    return Some($vertex_iter_type::new(vid, mesh));
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

            fn vertex_ids(&self) -> impl Iterator<Item = VertexId> {
                let mut he_iter = self.incoming_halfedges();
                let mut prev: Option<HalfedgeIter<'m, M>> = None;
                let mesh = unsafe { self.mesh.as_ref() };
                from_fn(move || {
                    loop {
                        if let Some(he_item) = prev.take() {
                            let hid = he_item.data.next();
                            let eid = mesh.he_edge(hid);
                            if mesh.e_halfedge(eid) == hid {
                                return Some(mesh.he_to(hid));
                            }
                        } else {
                            if let Some(he_item) = he_iter.next() {
                                let eid = mesh.he_edge(he_item.id);
                                if mesh.e_halfedge(eid) == he_item.id {
                                    let vid = mesh.he_to(he_item.data.prev());
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
    };
}

edge_vertex_iter_impl!(VertexIter, EdgeIter, VertexIter);
edge_vertex_iter_impl!(VertexIterMut, EdgeIter, VertexIter);

impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> VertexEdgesAndVerticesMut<'m, M>
    for VertexIterMut<'m, M>
{
    fn edges_mut(&mut self) -> impl Iterator<Item = EdgeIterMut<'m, M>> {
        let mut he_iter = self.incoming_halfedges_mut();
        let mut prev: Option<HalfedgeIterMut<'m, M>> = None;
        from_fn(move || {
            loop {
                if let Some(mut he_item) = prev.take() {
                    let mesh = unsafe { he_item.mesh.as_mut() };
                    let hid = he_item.data.next();
                    let eid = mesh.he_edge(hid);
                    if mesh.e_halfedge(eid) == hid {
                        return Some(EdgeIterMut::new(eid, mesh));
                    }
                } else {
                    if let Some(mut he_item) = he_iter.next() {
                        let mesh = unsafe { he_item.mesh.as_mut() };
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

    fn vertices_mut(&mut self) -> impl Iterator<Item = VertexIterMut<'m, M>> {
        let mut he_iter = self.incoming_halfedges_mut();
        let mut prev: Option<HalfedgeIterMut<'m, M>> = None;
        from_fn(move || {
            loop {
                if let Some(mut he_item) = prev.take() {
                    let mesh = unsafe { he_item.mesh.as_mut() };
                    let hid = he_item.data.next();
                    let eid = mesh.he_edge(hid);
                    if mesh.e_halfedge(eid) == hid {
                        return Some(VertexIterMut::new(mesh.he_to(hid), mesh));
                    }
                } else {
                    if let Some(mut he_item) = he_iter.next() {
                        let mesh = unsafe { he_item.mesh.as_mut() };
                        let eid = mesh.he_edge(he_item.id);
                        if mesh.e_halfedge(eid) == he_item.id {
                            let vid = mesh.he_to(he_item.data.prev());
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
