use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{
            Edge, EdgeHalfedge, EdgeMut, Halfedge, HalfedgeIterMethod, HalfedgeMut,
            HalfedgeIterMutMethod,
        },
        mesh::{Mesh, MeshCore},
    },
};

use super::{ElementId, HalfedgeData, HalfedgeId};

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

pub trait VertexData {
    fn halfedge(&self) -> HalfedgeId;
    fn set_halfedge(&mut self, hid: HalfedgeId);
}

element_iter_struct!(struct Vertex -> MeshCore, VertexId, VertexData, from, as_ref, vertex_data, {});
element_iter_struct!(struct VertexMut -> MeshCore, VertexId, VertexData, from_mut, as_mut, vertex_data_mut, { mut });

// Macro to generate circulate vertex structs
macro_rules! circulate_vertex_struct {
    ($name:ident, $halfedge_iter: ident, $halfedge_trait: ident, $halfedge_method: ident, $incoming_next:ident, {$( $mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            halfedge: $halfedge_iter<'m, M>,
            first: bool,
        }
        impl<'m, M: Mesh<VertexData: VertexData, HalfedgeData: HalfedgeData>> $name<'m, M> {

            #[inline]
            fn new(from_hid: HalfedgeId, mesh: &'m $($mut_)? M) -> Self {
                let hid = mesh.halfedge_data(from_hid).prev();
                Self {
                    first_hid: hid,
                    halfedge: $halfedge_iter::new(hid, mesh),
                    first: true,
                }
            }
        }

        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn valid(&self) -> bool {
                self.first || self.halfedge.id != self.first_hid
            }
        }

        impl <'m, M: Mesh<HalfedgeData : HalfedgeData>> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            fn next(&mut self) -> Option<Self::Item> where $halfedge_iter<'m, M>: $halfedge_trait<'m, M> {
                if !self.valid() {
                    return None;
                }
                self.first = false;
                let next_halfedge = self.halfedge.$incoming_next();
                Some(std::mem::replace(&mut self.halfedge, next_halfedge))
            }
        }
     };
}

circulate_vertex_struct!(
    CirculateVertex,
    Halfedge,
    HalfedgeIterMethod,
    halfedge,
    incoming_next,
    {}
);
circulate_vertex_struct!(CirculateVertexMut, HalfedgeMut, HalfedgeIterMutMethod, halfedge_mut, incoming_next_mut, { mut });

macro_rules! vertex_edges_struct {
    ($name:ident, $halfedge_iter: ident, $edge_iter: ident, $halfedge_trait: ident,  $incoming_next: ident, $halfedge_next: ident, $edge_method:ident, {$( $mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            halfedge: $halfedge_iter<'m, M>,
            visited: bool,
            first: bool,
        }

        impl <'m, M: Mesh<VertexData: VertexData, HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            fn new(vid: VertexId, mesh: &'m $($mut_)? M) -> Self {
                let next_hid = mesh.vertex_data(vid).halfedge();
                let hid = mesh.halfedge_data(next_hid).prev();
                Self {
                    first_hid: hid,
                    halfedge: $halfedge_iter::new(hid, mesh),
                    visited: false,
                    first: true,
                }
            }
        }

        impl <'m, M: Mesh<HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            fn next_halfedge(&mut self) {
                self.first = false;
                if self.visited {
                    self.visited = false;
                    self.halfedge = self.halfedge.$incoming_next();
                } else {
                    self.visited = true;
                }
            }
        }
        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn valid(&self) -> bool {
                self.first || self.visited || self.halfedge.id != self.first_hid
            }
        }

        impl <'m, M: Mesh<HalfedgeData : HalfedgeData>> Iterator for $name<'m, M> {
            type Item = $edge_iter<'m, M>;

            fn next(&mut self) -> Option<Self::Item> where $halfedge_iter<'m, M>: $halfedge_trait<'m, M> {
                loop {
                    if !self.valid() {
                        return None;
                    }
                    if self.visited {
                        let $($mut_)? halfedge = self.halfedge.$halfedge_next();
                        let edge = halfedge.$edge_method();
                        if edge.halfedge().id == halfedge.id {
                            self.next_halfedge();
                            return Some(edge);
                        }
                    } else {
                        let edge = self.halfedge.$edge_method();
                        if edge.halfedge().id == self.halfedge.id {
                            self.next_halfedge();
                            return Some(edge);
                        }
                    };
                    self.next_halfedge();
                }
            }
        }
    }
}
vertex_edges_struct!(
    VertexEdges,
    Halfedge,
    Edge,
    HalfedgeIterMethod,
    incoming_next,
    next,
    edge,
    {}
);
vertex_edges_struct!(VertexEdgesMut, HalfedgeMut, EdgeMut, HalfedgeIterMutMethod, incoming_next_mut, next_mut, edge_mut, {mut});

macro_rules! vertex_vertices_struct {
    ($name:ident, $halfedge_iter: ident, $edge_iter: ident, $vertex_iter: ident, $halfedge_trait: ident, $incoming_next: ident, $halfedge_next: ident, $edge_method: ident, $to_vertex: ident, $from_vertex: ident, {$( $mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            halfedge: $halfedge_iter<'m, M>,
            visited: bool,
            first: bool,
        }

        impl <'m, M: Mesh<VertexData: VertexData, HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            fn new(vid: VertexId, mesh: &'m $($mut_)? M) -> Self {
                let next_hid = mesh.vertex_data(vid).halfedge();
                let hid = mesh.halfedge_data(next_hid).prev();
                Self {
                    first_hid: hid,
                    halfedge: $halfedge_iter::new(hid, mesh),
                    visited: false,
                    first: true,
                }
            }
        }

        impl <'m, M: Mesh<HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            fn next_halfedge(&mut self) {
                self.first = false;
                if self.visited {
                    self.visited = false;
                    self.halfedge = self.halfedge.$incoming_next();
                } else {
                    self.visited = true;
                }
            }
        }

        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn valid(&self) -> bool {
                self.first || self.visited || self.halfedge.id != self.first_hid
            }
        }

        impl <'m, M: Mesh<HalfedgeData : HalfedgeData>> Iterator for $name<'m, M> {
            type Item = $vertex_iter<'m, M>;

            fn next(&mut self) -> Option<Self::Item> where $halfedge_iter<'m, M>: $halfedge_trait<'m, M> {
                loop {
                    if !self.valid() {
                        return None;
                    }

                    if self.visited {
                        let $($mut_)? halfedge = self.halfedge.$halfedge_next();
                        let edge = halfedge.$edge_method();
                        if edge.halfedge().id == halfedge.id {
                            self.next_halfedge();
                            return Some(halfedge.$to_vertex());
                        }
                    } else {
                        let edge = self.halfedge.edge();
                        if edge.halfedge().id == self.halfedge.id {
                            self.next_halfedge();
                            return Some(self.halfedge.$from_vertex());
                        }
                    };
                    self.next_halfedge();
                }
            }
        }

    }
}

vertex_vertices_struct!(
    VertexVertices,
    Halfedge,
    Edge,
    Vertex,
    HalfedgeIterMethod,
    incoming_next,
    next,
    edge,
    to,
    from,
    {}
);
vertex_vertices_struct!(VertexVerticesMut, HalfedgeMut, EdgeMut, VertexMut, HalfedgeIterMutMethod, incoming_next_mut, next_mut, edge_mut, to_mut, from_mut, {mut});

macro_rules! vertex_methods {
    ($name: ident, $incoming_method: ident, $circulate_vertex: ident, $vertex_edges: ident, $vertex_vertices: ident, $outgoing_method: ident, $halfedge_iter: ident, $halfedge_next: ident, $edges: ident, $vertices: ident, $into_ref: ident, {$( $mut_:tt )?}) => {
        impl<'m, M: Mesh<VertexData: VertexData, HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            pub fn $incoming_method(& $($mut_)? self) -> $circulate_vertex<'m, M> {
                unsafe { $circulate_vertex::new(self.data.halfedge(), self.mesh.$into_ref()) }
            }

            #[inline]
            pub fn $outgoing_method(& $($mut_)? self) -> impl Iterator<Item = $halfedge_iter<'m, M>> {
                self.$incoming_method().map(|$($mut_)? iter| iter.$halfedge_next())
            }

            #[inline]
            pub fn $edges(& $($mut_)? self) -> $vertex_edges<'m, M> {
                unsafe {$vertex_edges::new(self.id, self.mesh.$into_ref())}
            }

            #[inline]
            pub fn $vertices(& $($mut_)? self) -> $vertex_vertices<'m, M> {
                unsafe {$vertex_vertices::new(self.id, self.mesh.$into_ref())}
            }
        }
    };
}

vertex_methods!(
    Vertex,
    incoming_halfedges,
    CirculateVertex,
    VertexEdges,
    VertexVertices,
    outgoing_halfedges,
    Halfedge,
    next,
    edges,
    vertices,
    as_ref,
    {}
);
vertex_methods!(
    VertexMut,
    incoming_halfedges,
    CirculateVertex,
    VertexEdges,
    VertexVertices,
    outgoing_halfedges,
    Halfedge,
    next,
    edges,
    vertices,
    as_ref,
    {}
);
vertex_methods!(
    VertexMut,
    incoming_halfedges_mut,
    CirculateVertexMut,
    VertexEdgesMut,
    VertexVerticesMut,
    outgoing_halfedges_mut,
    HalfedgeMut,
    next_mut,
    edges_mut,
    vertices_mut,
    as_mut,
    {mut}
);
