use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{ HalfedgeData, HalfedgeId, Halfedge, HalfedgeMut, HalfedgeIterMethod, HalfedgeIterMutMethod},
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

pub trait EdgeData {
    fn halfedge(&self) -> HalfedgeId;
}

element_iter_struct!(struct Edge -> Mesh, EdgeId, EdgeData, from, as_ref, edge_data, {});
element_iter_struct!(struct EdgeMut -> Mesh, EdgeId, EdgeData, from_mut, as_mut, edge_data_mut, {mut});

pub trait EdgeHalfedge<'m, M: Mesh> {
    fn halfedge(&self) -> Halfedge<'m, M>;
}

pub trait EdgeHalfedgeMut<'m, M: Mesh>: EdgeHalfedge<'m, M> {
    fn halfedge_mut(&mut self) -> HalfedgeMut<'m, M>;
}

macro_rules! edge_halfedge_trait {
    ($name:ident, $halfedge_trait: ident, $halfedge_iter: ident, $halfedge_method: ident, $into_ref: ident, {$($mut_:tt )?}) => {
        impl <'m, M: Mesh> $halfedge_trait<'m, M> for $name<'m, M> {
            #[inline]
            default fn $halfedge_method(& $($mut_)? self) -> $halfedge_iter<'m, M> {
                unsafe {
                    let hid = self.mesh.as_ref().e_halfedge(self.id);
                    $halfedge_iter::new(hid, self.mesh.$into_ref())
                }
            }
        }

        impl <'m, M: Mesh<EdgeData: EdgeData>> $halfedge_trait<'m, M> for $name<'m, M> {
            #[inline]
            fn $halfedge_method(& $($mut_)? self) -> $halfedge_iter<'m, M> {
                unsafe {
                    $halfedge_iter::new(self.data.halfedge(), self.mesh.$into_ref())
                }
            }
        }
    }
}
edge_halfedge_trait!(Edge, EdgeHalfedge, Halfedge, halfedge, as_ref, {});
edge_halfedge_trait!(EdgeMut, EdgeHalfedge, Halfedge, halfedge, as_ref, {});
edge_halfedge_trait!(EdgeMut, EdgeHalfedgeMut, HalfedgeMut, halfedge_mut, as_mut, {mut});

macro_rules! edge_halfedges_struct {
    (struct $name:ident, $halfedge_method: ident, $halfedge_iter: ident, $halfedge_trait: ident, $sibling: ident, {$($mut_:tt )?}) => {
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

        impl<'m, M: Mesh<HalfedgeData: HalfedgeData>> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> where $halfedge_iter<'m, M>: $halfedge_trait<'m, M> {
                if !self.valid() {
                    return None;
                }
                self.first = false;
                let sibling = self.halfedge.$sibling();
                Some(std::mem::replace(&mut self.halfedge, sibling))
            }
        }
    };
}

edge_halfedges_struct!(struct EdgeHalfedges, halfedge, Halfedge, HalfedgeIterMethod, sibling, {});
edge_halfedges_struct!(struct EdgeHalfedgesMut, halfedge_mut, HalfedgeMut, HalfedgeIterMutMethod, sibling_mut, {mut});

macro_rules! edge_halfedges_method {
    (struct $name:ident, $halfedge_method: ident, $halfedges_method: ident, $halfedge_iter: ident, $edge_halfedges: ident, $into_ref:ident, {$($mut_:tt )?}) => {
        impl <'m, M: Mesh> $name<'m, M> {

            #[inline]
            pub fn $halfedges_method(& $($mut_)? self) -> $edge_halfedges<'m, M> {
                let halfedge = self.$halfedge_method();
                let first_hid = halfedge.id;
                $edge_halfedges {
                    first_hid,
                    halfedge,
                    first: true,
                }
            }
        }
    }
}

edge_halfedges_method!(struct Edge, halfedge, halfedges, Halfedge, EdgeHalfedges, as_ref, {});
edge_halfedges_method!(struct EdgeMut, halfedge, halfedges, Halfedge, EdgeHalfedges, as_ref, {});
edge_halfedges_method!(struct EdgeMut, halfedge_mut, halfedges_mut, HalfedgeMut, EdgeHalfedgesMut, as_mut, {mut});
