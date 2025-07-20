use std::{iter::from_fn, marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
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

element_iter_struct!(struct VertexIter -> MeshCore, VertexId, Vertex, from, as_ref, vertex, {});
element_iter_struct!(struct VertexIterMut -> MeshCore, VertexId, Vertex, from_mut, as_mut, vertex_mut, {mut});

trait HalfedgeIncomingNext {
    fn next(&mut self);
}

// Macro to generate circulate vertex structs
macro_rules! circulate_vertex_struct {
    ($name:ident, $from_ref:ident, $into_ref: ident, $he_data: ident, {$( $mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            hid: HalfedgeId,
            he: &'m $($mut_)? M::Halfedge,
            mesh: NonNull<M>,
            first: bool,
            _marker: PhantomData<&'m $($mut_)? M>,
        }
        impl<'m, M: Mesh<Vertex: Vertex, Halfedge: Halfedge>> $name<'m, M> {
            pub fn new(from_hid: HalfedgeId, mesh: &'m $($mut_)? M) -> Self {
                let $($mut_)? mesh = NonNull::$from_ref(mesh);
                unsafe {
                    let hid = mesh.as_ref().halfedge(from_hid).prev();
                    let he = mesh.$into_ref().$he_data(hid);
                    Self {
                        first: true,
                        first_hid: hid,
                        hid,
                        he,
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

        impl<'m, M: Mesh<Halfedge: Halfedge>> HalfedgeIncomingNext for $name<'m, M> {
            default fn next(&mut self) {
                unsafe {
                    self.hid = self.mesh.as_ref().he_incoming_next(self.hid);
                    self.he = self.mesh.$into_ref().$he_data(self.hid);
                    self.first = false;
                }
            }
        }

        impl<'m, M: Mesh<Halfedge: HalfedgeExt>> HalfedgeIncomingNext for $name<'m, M> {
            fn next(&mut self) {
                self.hid = self.he.incoming_next();
                unsafe {
                    self.he = self.mesh.$into_ref().$he_data(self.hid);
                }
                self.first = false;
            }
        }
    };
}

circulate_vertex_struct!(CirculateVertex, from, as_ref, halfedge, {});
circulate_vertex_struct!(CirculateVertexMut, from_mut, as_mut, halfedge_mut, {mut});

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
            pub fn halfedge(&self) -> HalfedgeIter<'m, M> {
                unsafe { HalfedgeIter::new(self.data.halfedge(), self.mesh.as_ref()) }
            }
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
