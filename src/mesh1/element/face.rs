use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    element_iter_struct, mesh1::{element::{HalfedgeIter, HalfedgeIterMut}, mesh::{Mesh, MeshCore}}, INVALID_IND
};

use super::{ElementId, Halfedge, HalfedgeId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

impl Default for FaceId {
    #[inline]
    fn default() -> Self {
        FaceId(INVALID_IND)
    }
}

impl From<usize> for FaceId {
    #[inline]
    fn from(index: usize) -> Self {
        FaceId(index)
    }
}

impl ElementId for FaceId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for FaceId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub trait Face {
    fn halfedge(&self) -> HalfedgeId;
    fn set_halfedge(&mut self, halfedge: HalfedgeId);
}

element_iter_struct!(struct FaceIter -> MeshCore, FaceId, Face, from, as_ref, face, {});
element_iter_struct!(struct FaceIterMut -> MeshCore, FaceId, Face, from_mut, as_mut, face_mut, {mut});


macro_rules! face_halfedges_struct {
    ($name:ident, $halfedge_iter: ident, $halfedge_next: ident, $halfedge_prev: ident,  $into_ref:ident, {$($mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            halfedge: $halfedge_iter<'m, M>,
            first: bool,
        }

        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn new(hid: HalfedgeId, mesh: &'m $($mut_)? M) -> Self {
                Self {
                    first_hid: hid,
                    halfedge: $halfedge_iter::new(hid, mesh),
                    first: true,
                }
            }

            #[inline]
            fn valid(&self) -> bool {
                self.first || self.halfedge.id != self.first_hid
            }
        }

        impl<'m, M: Mesh<Halfedge: Halfedge>> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }

                self.first = false;
                let he_next = self.halfedge.$halfedge_next();
                Some(std::mem::replace(&mut self.halfedge, he_next))
            }
        }

        impl<'m, M: Mesh<Halfedge: Halfedge>> DoubleEndedIterator for $name<'m, M> {
            #[inline]
            fn next_back(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }

                self.first = false;
                let he_prev = self.halfedge.$halfedge_prev();
                Some(std::mem::replace(&mut self.halfedge, he_prev))
            }
        }
    }
}

face_halfedges_struct!(FaceHalfedges, HalfedgeIter, next, prev,  as_ref, {});
face_halfedges_struct!(FaceHalfedgesMut, HalfedgeIterMut, next_mut, prev_mut,  as_mut, {mut});

macro_rules! face_halfedges_method {
    (struct $name:ident, $halfedge_method: ident, $halfedges_method: ident, $halfedge_iter: ident, $edge_halfedges: ident, $into_ref:ident, {$($mut_:tt )?}) => {
        impl <'m, M: Mesh<Halfedge: Halfedge, Face: Face>> $name<'m, M> {
            #[inline]
            pub fn $halfedge_method(& $($mut_)? self) -> $halfedge_iter<'m, M> {
                unsafe {
                    $halfedge_iter::new(self.data.halfedge(), self.mesh.$into_ref())
                }
            }

            #[inline]
            pub fn $halfedges_method(& $($mut_)? self) -> $edge_halfedges<'m, M> {
                unsafe {
                    $edge_halfedges::new(self.data.halfedge(), self.mesh.$into_ref())
                }
            }
        }
    }
}

face_halfedges_method!(struct FaceIter, halfedge, halfedges, HalfedgeIter, FaceHalfedges, as_ref, {});
face_halfedges_method!(struct FaceIterMut, halfedge, halfedges, HalfedgeIter, FaceHalfedges, as_ref, {});
face_halfedges_method!(struct FaceIterMut, halfedge_mut, halfedges_mut, HalfedgeIterMut, FaceHalfedgesMut, as_mut, {mut});
