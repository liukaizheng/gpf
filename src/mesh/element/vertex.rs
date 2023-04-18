use std::ops::{Deref, Index};

use super::{iter_next, Edge, Element, Halfedge, HalfedgeId, HalfedgeIter};
use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(usize);

impl From<usize> for VertexId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl<T> Index<VertexId> for Vec<T> {
    type Output = T;
    fn index(&self, index: VertexId) -> &Self::Output {
        &self[index.0]
    }
}

pub struct VertexIter<'a, M: Mesh> {
    id: VertexId,
    mesh: &'a M,
}

impl<'a, M: Mesh> VertexIter<'a, M> {
    pub fn new(id: VertexId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

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
    fn halfedge(&self) -> HalfedgeIter<'_, Self::M>;
}

impl<'a, M: Mesh> Vertex for VertexIter<'a, M> {
    fn halfedge(&self) -> HalfedgeIter<'_, Self::M> {
        HalfedgeIter::new(self.mesh.v_halfedge(self.id), self.mesh)
    }
}

impl<'a, M: Mesh> Element for VertexIter<'a, M> {
    type Item = VertexId;
    type M = M;

    fn id(&self) -> VertexId {
        self.id
    }
    fn mesh(&self) -> &M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.vertex_is_valid(self.id)
    }

    fn next(&mut self) -> bool {
        self.id.0 += 1;
        self.id.0 < self.capacity()
    }
}

impl<'a, M: Mesh> Iterator for VertexIter<'a, M> {
    type Item = VertexId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
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
	curr_state: VNeighborIterState,
	start_state: VNeighborIterState,
	first_he: HalfedgeId
}

impl <'a, M: Mesh> VVertexIter<'a, M> {
	fn new(id: VertexId, mesh: &'a M) -> Self {
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
		Self { mesh, curr_state, start_state, first_he }
	}
}

impl <'a, M: Mesh> Element for VVertexIter<'a, M> {
	type Item = VertexId;
	type M = M;

	#[inline(always)]
	fn id(&self) -> VertexId {
		if self.curr_state.processing_incoming {
			self.mesh.he_vertex(self.curr_state.curr_he)
		} else {
			self.mesh.he_tip_vertex(self.curr_state.curr_he)
		}
	}

	fn mesh(&self) -> &M {
		self.mesh
	}

	fn valid(&self) -> bool {
		self.curr_state.is_halfedge_canonical(self.mesh)
	}

	fn next(&mut self) -> bool {
		self.curr_state.next(self.mesh, &mut self.first_he);
		self.curr_state != self.start_state
	}
}

impl <'a, M: Mesh> Iterator for VVertexIter<'a, M> {
	type Item = VertexId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}