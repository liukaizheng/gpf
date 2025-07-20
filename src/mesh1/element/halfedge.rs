use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{
            EdgeId, EdgeIter, EdgeIterMut, FaceIter, FaceIterMut, VertexId, VertexIter,
            VertexIterMut,
        },
        mesh::{Mesh, MeshCore},
    },
};

use super::{ElementId, FaceId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

pub trait Halfedge {
    fn vertex(&self) -> VertexId;
    fn set_vertex(&mut self, vertex: VertexId);

    fn prev(&self) -> HalfedgeId;
    fn set_prev(&mut self, prev: HalfedgeId);

    fn next(&self) -> HalfedgeId;
    fn set_next(&mut self, next: HalfedgeId);

    fn face(&self) -> FaceId;
    fn set_face(&mut self, face: FaceId);
}

pub trait HalfedgeExt: Halfedge {
    fn edge(&self) -> EdgeId;
    fn sibling(&self) -> HalfedgeId;
    fn incoming_next(&self) -> HalfedgeId;
}

element_iter_struct!(struct HalfedgeIter -> MeshCore, HalfedgeId, Halfedge, from, as_ref, halfedge, {});
element_iter_struct!(struct HalfedgeIterMut -> MeshCore, HalfedgeId, Halfedge, from_mut, as_mut, halfedge_mut, {mut});

macro_rules! halfedge_base_methods {
    (struct $name:ident, $vertex: ident, $halfedge: ident, $face: ident, $from: ident, $to: ident, $next: ident, $prev: ident, $face_method: ident, $into_ref:ident, {$( $mut_:tt )?}) => {
        impl<'m, M: MeshCore<Halfedge: Halfedge>> $name<'m, M> {
            #[inline]
            fn $from(& $($mut_)? self) -> $vertex<'m, M> {
                unsafe {
                    let vid = self.mesh.as_ref().he_to(self.data.prev());
                    $vertex::new(vid, self.mesh.$into_ref())
                }
            }
            #[inline]
            fn $to(& $($mut_)? self) -> $vertex<'m, M> {
                unsafe { $vertex::new(self.data.vertex(), self.mesh.$into_ref()) }
            }

            #[inline]
            fn $next(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe { $halfedge::new(self.data.next(), self.mesh.$into_ref()) }
            }

            #[inline]
            fn $prev(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe { $halfedge::new(self.data.prev(), self.mesh.$into_ref()) }
            }

            #[inline]
            fn $face_method(& $($mut_)? self) -> $face<'m, M> {
                unsafe { $face::new(self.data.face(), self.mesh.$into_ref()) }
            }
        }
    };
}
halfedge_base_methods! (struct HalfedgeIter, VertexIter, HalfedgeIter, FaceIter, from, to, next, prev, face, as_ref, {});
halfedge_base_methods! (struct HalfedgeIterMut, VertexIter, HalfedgeIter, FaceIter, from, to, next, prev, face, as_ref, {});
halfedge_base_methods! (struct HalfedgeIterMut, VertexIterMut, HalfedgeIterMut, FaceIterMut, from_mut, to_mut, next_mut, prev_mut, face_mut, as_mut, { mut });

pub trait HalfedgeIterMethod<'m, M: Mesh> {
    fn edge(&self) -> EdgeIter<'m, M>;
    fn sibling(&self) -> HalfedgeIter<'m, M>;
    fn incoming_next(&self) -> HalfedgeIter<'m, M>;
}

pub trait HalfedgeIterMutMethod<'m, M: Mesh>: HalfedgeIterMethod<'m, M> {
    fn edge_mut(&mut self) -> EdgeIterMut<'m, M>;
    fn sibling_mut(&mut self) -> HalfedgeIterMut<'m, M>;
    fn incoming_next_mut(&mut self) -> HalfedgeIterMut<'m, M>;
}

macro_rules! halfedge_base_methods {
    (struct $name:ident -> $halfedge_trait: ident, $halfedge: ident, $edge: ident, $edge_method: ident, $sibling: ident, $incoming_next: ident, $into_ref:ident, {$( $mut_:tt )?}) => {
        default impl<'m, M: Mesh<Halfedge: Halfedge>> $halfedge_trait<'m, M> for $name<'m, M> {
            #[inline]
            fn $edge_method(& $($mut_)? self) -> $edge<'m, M> {
                unsafe {
                    let eid = self.mesh.as_ref().he_edge(self.id);
                    $edge::new(eid, self.mesh.$into_ref())
                }
            }

            #[inline]
            fn $sibling(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.mesh.as_ref().he_sibling(self.id);
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }

            #[inline]
            fn $incoming_next(& $($mut_)? self) -> $halfedge<'m, M> {
                unsafe {
                    let hid = self.mesh.as_ref().he_incoming_next(self.id);
                    $halfedge::new(hid, self.mesh.$into_ref())
                }
            }
        }

        impl<'m, M: Mesh<Halfedge: HalfedgeExt>> $halfedge_trait<'m, M> for $name<'m, M> {
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

halfedge_base_methods! (struct HalfedgeIter -> HalfedgeIterMethod, HalfedgeIter, EdgeIter, edge, sibling, incoming_next, as_ref, {});
halfedge_base_methods! (struct HalfedgeIterMut -> HalfedgeIterMethod, HalfedgeIter, EdgeIter, edge, sibling, incoming_next, as_ref, {});
halfedge_base_methods! (struct HalfedgeIterMut -> HalfedgeIterMutMethod, HalfedgeIterMut, EdgeIterMut, edge_mut, sibling_mut, incoming_next_mut, as_mut, {mut});
