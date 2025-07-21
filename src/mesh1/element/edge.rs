use std::{marker::PhantomData, ops::Deref, ptr::NonNull};

use crate::{
    INVALID_IND, element_iter_struct,
    mesh1::{
        element::{Halfedge, HalfedgeExt, HalfedgeId, HalfedgeIter},
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

macro_rules! edge_halfedges_struct {
    (struct $name:ident, $halfedge: ident, $halfedge_method: ident, $halfedge_iter: ident, $from_ref:ident, $into_ref:ident, {$($mut_:tt )?}, { $( $set_hid:tt )* }) => {
        struct $name<'m, M: Mesh> {
            first_hid: HalfedgeId,
            hid: HalfedgeId,
            he: &'m $($mut_)? M::$halfedge,
            mesh: NonNull<M>,
            first: bool,
            _marker: PhantomData<&'m M>,
        }

        impl<'m, M: Mesh> $name<'m, M> {
            #[inline]
            fn new(first_hid: HalfedgeId, mesh: &'m M) -> Self {
                unsafe {
                    let mesh = NonNull::$from_ref(mesh);
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

        impl<'m, M: Mesh> Iterator for $name<'m, M> {
            type Item = $halfedge_iter<'m, M>;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                if !self.valid() {
                    return None;
                }
                unsafe {
                    self.first = false;
                    let ret = $halfedge_iter::new_with_data(self.hid, self.he, self.mesh.$into_ref());
                    $($set_hid)*
                    Some(ret)
                }
            }
        }
    };
}

edge_halfedges_struct!(struct EdgeHalfedges, Halfedge, halfedge, HalfedgeIter, from, as_ref, {}, {
    (self.hid = self.mesh.as_ref().he_sibling(self.hid));
});

struct EdgeHalfedgesFast<'m, M: Mesh> {
    first_hid: HalfedgeId,
    hid: HalfedgeId,
    he: &'m M::Halfedge,
    mesh: &'m M,
    first: bool,
    _marker: PhantomData<&'m M>,
}

impl<'m, M: Mesh> EdgeHalfedgesFast<'m, M> {
    #[inline]
    fn new(first_hid: HalfedgeId, mesh: &'m M) -> Self {
        Self {
            first_hid,
            hid: first_hid,
            he: mesh.halfedge(first_hid),
            mesh,
            first: true,
            _marker: PhantomData,
        }
    }

    #[inline]
    fn valid(&self) -> bool {
        self.first || self.hid != self.first_hid
    }
}

impl<'m, M: Mesh<Halfedge: HalfedgeExt>> Iterator for EdgeHalfedgesFast<'m, M> {
    type Item = HalfedgeIter<'m, M>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if !self.valid() {
            return None;
        }
        self.first = false;
        let ret = HalfedgeIter::new(self.hid, self.mesh);
        self.hid = self.he.sibling();
        Some(ret)
    }
}

pub trait EdgeMethod<'m, M: Mesh> {
    fn halfedges(&self) -> EdgeHalfedges<'m, M>;
}

pub trait EdgeMethodFast<'m, M: Mesh<Halfedge: HalfedgeExt>> {
    fn halfedges(&self) -> EdgeHalfedgesFast<'m, M>;
}

impl<'m, M: Mesh> EdgeMethod<'m, M> for EdgeIter<'m, M> {
    fn halfedges(&self) -> EdgeHalfedges<'m, M> {
        unsafe { EdgeHalfedges::new(self.mesh.as_ref().e_halfedge(self.id), self.mesh.as_ref()) }
    }
}

impl<'m, M: Mesh<Edge: Edge, Halfedge: HalfedgeExt>> EdgeMethodFast<'m, M> for EdgeIter<'m, M> {
    fn halfedges(&self) -> EdgeHalfedgesFast<'m, M> {
        unsafe { EdgeHalfedgesFast::new(self.data.halfedge(), self.mesh.as_ref()) }
    }
}
