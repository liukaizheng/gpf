use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Edge, EdgeId, Element, ElementId, Halfedge, HalfedgeId, HalfedgeIter};
use crate::{element_id, element_iterator, halfedges_iterator, mesh::Mesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(pub usize);

element_id! {struct VertexId}

pub struct VertexIter<'a, M: Mesh> {
    id: VertexId,
    mesh: &'a M,
}

impl<'a, M: Mesh> VertexIter<'a, M> {
    pub fn new(id: VertexId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_vertices()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_vertices_capacity()
    }
}

impl<'a, M: Mesh> Deref for VertexIter<'a, M> {
    type Target = VertexId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

pub trait Vertex: Element {
    fn halfedge(&self) -> HalfedgeIter<Self::M>;
    fn adjacent_vertices(&self) -> VVertexIter<Self::M>;
    fn adjacent_edges(&self) -> VEdgeIter<Self::M>;
    fn incoming_halfedge(&self) -> VIncomingHalfedgeIter<Self::M>;
    fn outgoing_halfedge(&self) -> VOutgoingHalfedgeIter<Self::M>;
}

impl<'a, M: Mesh> Vertex for VertexIter<'a, M> {
    #[inline(always)]
    fn halfedge(&self) -> HalfedgeIter<Self::M> {
        HalfedgeIter::new(self.mesh.v_halfedge(self.id), self.mesh)
    }

    #[inline(always)]
    fn adjacent_vertices(&self) -> VVertexIter<Self::M> {
        VVertexIter::new(self.id(), self.mesh)
    }

    #[inline(always)]
    fn adjacent_edges(&self) -> VEdgeIter<Self::M> {
        VEdgeIter::new(self.id, self.mesh)
    }

    #[inline(always)]
    fn incoming_halfedge(&self) -> VIncomingHalfedgeIter<Self::M> {
        let mut he = self.halfedge();
        he.prev();
        VIncomingHalfedgeIter::new(self.mesh, *he)
    }

    #[inline(always)]
    fn outgoing_halfedge(&self) -> VOutgoingHalfedgeIter<Self::M> {
        VOutgoingHalfedgeIter::new(self.mesh, self.mesh.v_halfedge(self.id))
    }
}

element_iterator! {
    struct VertexIter -> VertexId, {
        fn id(&self) -> VertexId {
            self.id
        }
    }, {
        fn valid(&self) -> bool {
            self.mesh.vertex_is_valid(self.id)
        }
    }, {
        fn next(&mut self) {
            *self.id += 1;
        }
    }, {
        fn is_end(&self) -> bool {
            *self.id == self.capacity()
        }
    }
}

/// adjacent vertices
#[derive(Clone, PartialEq, Eq)]
pub struct VNeighborIterState {
    processing_incoming: bool,
    curr_he: HalfedgeId,
}

impl VNeighborIterState {
    pub fn new(he: HalfedgeId) -> Self {
        Self {
            processing_incoming: false,
            curr_he: he,
        }
    }

    #[inline]
    pub fn is_halfedge_canonical<M: Mesh>(&self, mesh: &M) -> bool {
        if mesh.use_implicit_twin() {
            true
        } else {
            self.curr_he == *mesh.halfedge(self.curr_he).edge().halfedge()
        }
    }

    #[inline]
    pub fn next<M: Mesh>(&mut self, mesh: &M, first_he: &mut HalfedgeId) {
        if mesh.use_implicit_twin() {
            self.curr_he = mesh.he_next_outgoing_neighbor(self.curr_he);
        } else {
            if self.processing_incoming {
                self.curr_he = mesh.he_next_incoming_neighbor(self.curr_he);
                if self.curr_he == *first_he {
                    self.processing_incoming = false;
                    self.curr_he = mesh.he_next(self.curr_he);
                    *first_he = self.curr_he;
                }
            } else {
                self.curr_he = mesh.he_next_outgoing_neighbor(self.curr_he);
                if self.curr_he == *first_he {
                    self.processing_incoming = true;
                    self.curr_he = mesh.he_prev(self.curr_he);
                    *first_he = self.curr_he;
                }
            }
        }
    }
}

pub struct VVertexIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    curr_state: VNeighborIterState,
    start_state: VNeighborIterState,
    first_he: HalfedgeId,
}

macro_rules! vertex_state_iterator {
    (struct $name:ident -> $item: ty, {$($id: tt)*}) => {
        impl<'a, M: Mesh> $name<'a, M> {
            pub fn new(id: VertexId, mesh: &'a M) -> Self {
                let mut first_he = mesh.v_halfedge(id);
                let mut curr_state = VNeighborIterState::new(first_he);
                let first_state = curr_state.clone();
                while !curr_state.is_halfedge_canonical(mesh) {
                    curr_state.next(mesh, &mut first_he);
                    if curr_state == first_state {
                        break;
                    }
                }
                let start_state = curr_state.clone();
                Self {
                    mesh,
                    just_start: true,
                    curr_state,
                    start_state,
                    first_he,
                }
            }
        }
        element_iterator!{struct $name -> $item, {
            $($id)*
        }, {
            fn valid(&self) -> bool {
                self.curr_state.is_halfedge_canonical(self.mesh)
            }
        }, {
            fn next(&mut self) {
                self.just_start = false;
                self.curr_state.next(self.mesh, &mut self.first_he);
            }
        }, {
            fn is_end(&self) -> bool {
                !self.just_start && self.curr_state == self.start_state
            }
        }}
    };
}

vertex_state_iterator! { struct VVertexIter -> VertexId, {
    fn id(&self) -> VertexId {
        if self.curr_state.processing_incoming {
            self.mesh.he_vertex(self.curr_state.curr_he)
        } else {
            self.mesh.he_tip_vertex(self.curr_state.curr_he)
        }
    }
}}

pub struct VIncomingHalfedgeIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, M: Mesh> VIncomingHalfedgeIter<'a, M> {
    #[inline(always)]
    pub fn new(mesh: &'a M, he: HalfedgeId) -> Self {
        Self {
            mesh,
            just_start: true,
            first_he: he,
            curr_he: he,
        }
    }

    #[inline(always)]
    fn next(&mut self) {
        self.curr_he = self.mesh.he_next_incoming_neighbor(self.curr_he);
    }
}

halfedges_iterator! {struct VIncomingHalfedgeIter}

pub struct VOutgoingHalfedgeIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    first_he: HalfedgeId,
    curr_he: HalfedgeId,
}

impl<'a, M: Mesh> VOutgoingHalfedgeIter<'a, M> {
    #[inline(always)]
    pub fn new(mesh: &'a M, he: HalfedgeId) -> Self {
        Self {
            mesh,
            just_start: true,
            first_he: he,
            curr_he: he,
        }
    }

    #[inline(always)]
    fn next(&mut self) {
        self.curr_he = self.mesh.he_next_outgoing_neighbor(self.curr_he);
    }
}

halfedges_iterator! {struct VOutgoingHalfedgeIter}

pub struct VEdgeIter<'a, M: Mesh> {
    mesh: &'a M,
    just_start: bool,
    curr_state: VNeighborIterState,
    start_state: VNeighborIterState,
    first_he: HalfedgeId,
}

vertex_state_iterator! { struct VEdgeIter -> EdgeId, {
    fn id(&self) -> EdgeId {
        self.mesh.he_edge(self.curr_state.curr_he)
    }
}}
