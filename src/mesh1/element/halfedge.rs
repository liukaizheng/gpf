use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND,
    mesh1::{
        element::{Edge, EdgeId, EdgeMut, Face, FaceMut, Vertex, VertexId, VertexMut},
        mesh::{Mesh, MeshCore},
    },
};

use super::{ElementId, FaceId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HalfedgeId(pub usize);

impl Default for HalfedgeId {
    #[inline]
    fn default() -> Self {
        HalfedgeId(INVALID_IND)
    }
}

impl From<usize> for HalfedgeId {
    #[inline]
    fn from(index: usize) -> Self {
        HalfedgeId(index)
    }
}

impl ElementId for HalfedgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for HalfedgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait HalfedgeData {
    fn vertex(&self) -> VertexId;
    fn set_vertex(&mut self, vertex: VertexId);

    fn prev(&self) -> HalfedgeId;
    fn set_prev(&mut self, prev: HalfedgeId);

    fn next(&self) -> HalfedgeId;
    fn set_next(&mut self, next: HalfedgeId);

    fn face(&self) -> FaceId;
    fn set_face(&mut self, face: FaceId);
}

pub trait HalfedgeDataExt: HalfedgeData {
    fn edge(&self) -> EdgeId;
    fn sibling(&self) -> HalfedgeId;
    fn incoming_next(&self) -> HalfedgeId;
}

element_handle_struct!(struct Halfedge -> MeshCore, HalfedgeId, HalfedgeData, from, as_ref, halfedge_data, {});
element_handle_struct!(struct HalfedgeMut -> MeshCore, HalfedgeId, HalfedgeData, from_mut, as_mut, halfedge_data_mut, {mut});

macro_rules! impl_halfedge_base_methods {
    (struct $name:ident, $vertex: ident, $halfedge: ident, $face: ident, $from: ident, $to: ident, $next: ident, $prev: ident, $face_method: ident, $into_ref:ident, {$( $mut_:tt )?}) => {
        impl<'m, M: MeshCore<HalfedgeData: HalfedgeData>> $name<'m, M> {
            #[inline]
            pub fn $from(& $($mut_)? self) -> $vertex<'m, M> {
                unsafe {
                    let vid = self.mesh.as_ref().he_to(self.data.prev());
                    $vertex::new(vid, self.mesh.$into_ref())
                }
            }
            #[inline]
            pub fn $to(& $($mut_)? self) -> $vertex<'m, M> {
                unsafe { $vertex::new(self.data.vertex(), self.mesh.$into_ref()) }
            }

            #[inline]
            pub fn $next(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe { $halfedge::new(self.data.next(), self.mesh.$into_ref()) }
            }

            #[inline]
            pub fn $prev(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe { $halfedge::new(self.data.prev(), self.mesh.$into_ref()) }
            }

            #[inline]
            pub fn $face_method(& $($mut_)? self) -> $face<'m, M> {
                unsafe { $face::new(self.data.face(), self.mesh.$into_ref()) }
            }
        }
    };
}
impl_halfedge_base_methods! (struct Halfedge, Vertex, Halfedge, Face, from, to, next, prev, face, as_ref, {});
impl_halfedge_base_methods! (struct HalfedgeMut, Vertex, Halfedge, Face, from, to, next, prev, face, as_ref, {});
impl_halfedge_base_methods! (struct HalfedgeMut, VertexMut, HalfedgeMut, FaceMut, from_mut, to_mut, next_mut, prev_mut, face_mut, as_mut, { mut });

pub trait HalfedgeNavigation<'m, M: Mesh> {
    fn edge(&self) -> Edge<'m, M>;
    fn sibling(&self) -> Halfedge<'m, M>;
    fn incoming_next(&self) -> Halfedge<'m, M>;
}

pub trait HalfedgeNavigationMut<'m, M: Mesh>: HalfedgeNavigation<'m, M> {
    fn edge_mut(&mut self) -> EdgeMut<'m, M>;
    fn sibling_mut(&mut self) -> HalfedgeMut<'m, M>;
    fn incoming_next_mut(&mut self) -> HalfedgeMut<'m, M>;
}

macro_rules! impl_halfedge_base_methods {
    (struct $name:ident -> $halfedge_trait: ident, $halfedge: ident, $edge: ident, $edge_method: ident, $sibling: ident, $incoming_next: ident, $into_ref:ident, {$( $mut_:tt )?}) => {
        impl<'m, M: Mesh<HalfedgeData: HalfedgeData>> $halfedge_trait<'m, M> for $name<'m, M> {
            #[inline]
            default fn $edge_method(& $($mut_)? self) -> $edge<'m, M> {
                unsafe {
                    let eid = self.mesh.as_ref().he_edge(self.id);
                    $edge::new(eid, self.mesh.$into_ref())
                }
            }

            #[inline]
            default fn $sibling(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.mesh.as_ref().he_sibling(self.id);
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }

            #[inline]
            default fn $incoming_next(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.mesh.as_ref().he_incoming_next(self.id);
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }
        }

        impl<'m, M: Mesh<HalfedgeData: HalfedgeDataExt>> $halfedge_trait<'m, M> for $name<'m, M> {
            #[inline]
            fn $edge_method(& $($mut_)? self) -> $edge<'m, M> {
                unsafe {
                    let eid = self.data.edge();
                    $edge::new(eid, self.mesh.$into_ref())
                }
            }

            #[inline]
            fn $sibling(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.data.sibling();
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }

            #[inline]
            fn $incoming_next(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.data.incoming_next();
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }
        }
    };
}

impl_halfedge_base_methods! (struct Halfedge -> HalfedgeNavigation, Halfedge, Edge, edge, sibling, incoming_next, as_ref, {});
impl_halfedge_base_methods! (struct HalfedgeMut -> HalfedgeNavigation, Halfedge, Edge, edge, sibling, incoming_next, as_ref, {});
impl_halfedge_base_methods! (struct HalfedgeMut -> HalfedgeNavigationMut, HalfedgeMut, EdgeMut, edge_mut, sibling_mut, incoming_next_mut, as_mut, {mut});

impl<'m, M: Mesh<HalfedgeData: HalfedgeData>> HalfedgeMut<'m, M> {
    #[inline]
    pub fn connect(&mut self, other: &mut Self) {
        self.data.set_next(other.id);
        other.data.set_prev(self.id);
    }
}
