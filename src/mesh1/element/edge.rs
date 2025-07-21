use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{ HalfedgeExt, HalfedgeId, HalfedgeIter, HalfedgeIterMut},
        mesh::Mesh,
    },
};

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

impl Default for EdgeId {
    #[inline]
    fn default() -> Self {
        EdgeId(INVALID_IND)
    }
}

impl From<usize> for EdgeId {
    #[inline]
    fn from(index: usize) -> Self {
        EdgeId(index)
    }
}

impl ElementId for EdgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for EdgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Edge {
    fn halfedge(&self) -> HalfedgeId;
}

element_iter_struct!(struct EdgeIter -> Mesh, EdgeId, Edge, from, as_ref, edge, {});
element_iter_struct!(struct EdgeIterMut -> Mesh, EdgeId, Edge, from_mut, as_mut, edge_mut, {mut});

trait EdgeHalfedge<'m, M: Mesh> {
    fn one_halfedge(&self) -> HalfedgeId;
}

macro_rules! edge_halfedge_trait {
    ($name:ident) => {
        impl<'m, M: Mesh> EdgeHalfedge<'m, M> for $name<'m, M> {
            #[inline]
            default fn one_halfedge(&self) -> HalfedgeId {
                unsafe { self.mesh.as_ref().e_halfedge(self.id) }
            }
        }

        impl<'m, M: Mesh<Edge: Edge>> EdgeHalfedge<'m, M> for $name<'m, M> {
            #[inline]
            fn one_halfedge(&self) -> HalfedgeId {
                self.data.halfedge()
            }
        }
    };
}

edge_halfedge_trait!(EdgeIter);
edge_halfedge_trait!(EdgeIterMut);

trait HalfedgeOperation<'m, M: Mesh> {
    fn next_halfedge(&mut self);
}

macro_rules! edge_halfedges_struct {
    (struct $name:ident, $halfedge: ident, $halfedge_method: ident, $halfedge_iter: ident, $from_ref:ident, $into_ref:ident, {$($mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            hid: HalfedgeId,
            he: &'m $($mut_)? M::$halfedge,
            mesh: NonNull<M>,
            first: bool,
            _marker: PhantomData<&'m M>,
        }

        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn new(first_hid: HalfedgeId, mesh: &'m $($mut_)? M) -> Self {
                unsafe {
                    let $($mut_)? mesh = NonNull::$from_ref(mesh);
                    let he = mesh.$into_ref().$halfedge_method(first_hid);
                    Self {
                        first_hid,
                        hid: first_hid,
                        he,
                        mesh,
                        first: true,
                        _marker: PhantomData,
                    }
                }
            }

            #[inline]
            fn valid(&self) -> bool {
                self.first || self.hid != self.first_hid
            }
        }

        impl <'m, M: Mesh> HalfedgeOperation<'m, M> for $name<'m, M> {
            #[inline]
            default fn next_halfedge(&mut self) {
                unsafe {
                    self.hid = self.mesh.as_ref().he_sibling(self.hid);
                }
            }
        }

        impl <'m, M: Mesh<Halfedge : HalfedgeExt>> HalfedgeOperation<'m, M> for $name<'m, M> {
            #[inline]
            fn next_halfedge(&mut self) {
                self.hid = self.he.sibling();
            }
        }

        impl<'m, M: Mesh> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }
                unsafe {
                    self.first = false;
                    let old_hid = self.hid;
                    self.next_halfedge();
                    let he = std::mem::replace(&mut self.he, self.mesh.$into_ref().$halfedge_method(self.hid));
                    let ret = $halfedge_iter::new_with_data(old_hid, he, self.mesh.$into_ref());
                    Some(ret)
                }
            }
        }
    };
}

edge_halfedges_struct!(struct EdgeHalfedges, Halfedge, halfedge, HalfedgeIter, from, as_ref, {});
edge_halfedges_struct!(struct EdgeHalfedgesMut, Halfedge, halfedge_mut, HalfedgeIterMut, from_mut, as_mut, {mut});

macro_rules! edge_halfedges_method {
    (struct $name:ident, $halfedge_method: ident, $halfedges_method: ident, $halfedge_iter: ident, $edge_halfedges: ident, $into_ref:ident, {$($mut_:tt )?}) => {
        impl <'m, M: Mesh> $name<'m, M> {
            #[inline]
            pub fn $halfedge_method(& $($mut_)? self) -> $halfedge_iter<'m, M> {
                unsafe {
                    $halfedge_iter::new(self.one_halfedge(), self.mesh.$into_ref())
                }
            }

            #[inline]
            pub fn $halfedges_method(& $($mut_)? self) -> $edge_halfedges<'m, M> {
                unsafe {
                    $edge_halfedges::new(self.one_halfedge(), self.mesh.$into_ref())
                }
            }
        }
    }
}

edge_halfedges_method!(struct EdgeIter, halfedge, halfedges, HalfedgeIter, EdgeHalfedges, as_ref, {});
edge_halfedges_method!(struct EdgeIterMut, halfedge, halfedges, HalfedgeIter, EdgeHalfedges, as_ref, {});
edge_halfedges_method!(struct EdgeIterMut, halfedge_mut, halfedges_mut, HalfedgeIterMut, EdgeHalfedgesMut, as_mut, {mut});
