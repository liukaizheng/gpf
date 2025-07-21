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
    (struct $name:ident, $halfedge_method: ident, $halfedge_iter: ident, $from_ref:ident, $into_ref:ident, {$($mut_:tt )?}) => {
        pub struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            hid: HalfedgeId,
            he: &'m $($mut_)? M::Halfedge,
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

        impl<'m, M: Mesh<Halfedge: Halfedge>> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }
                unsafe {
                    self.first = false;
                    let old_hid = self.hid;
                    self.hid = self.he.next();
                    let he = std::mem::replace(&mut self.he, self.mesh.$into_ref().$halfedge_method(self.hid));
                    let ret = $halfedge_iter::new_with_data(old_hid, he, self.mesh.$into_ref());
                    Some(ret)
                }
            }
        }

        impl<'m, M: Mesh<Halfedge: Halfedge>> DoubleEndedIterator for $name<'m, M> {
            fn next_back(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }
                unsafe {
                    self.first = false;
                    let old_hid = self.hid;
                    self.hid = self.he.prev();
                    let he = std::mem::replace(&mut self.he, self.mesh.$into_ref().$halfedge_method(self.hid));
                    let ret = $halfedge_iter::new_with_data(old_hid, he, self.mesh.$into_ref());
                    Some(ret)
                }
            }
        }
    }
}

face_halfedges_struct!(struct FaceHalfedges, halfedge, HalfedgeIter, from, as_ref, {});
face_halfedges_struct!(struct FaceHalfedgesMut, halfedge_mut, HalfedgeIterMut, from_mut, as_mut, {mut});

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
