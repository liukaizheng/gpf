use std::ops::{Index, Deref};

use crate::mesh::Mesh;
use super::{Element, iter_next, EdgeIter};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(usize);

impl From<usize> for HalfedgeId {
	fn from(id: usize) -> Self {
		Self(id)
	}
}

impl <T> Index<HalfedgeId> for Vec<T> {
	type Output = T;
	fn index(&self, index: HalfedgeId) -> &Self::Output {
		&self[index.0]
	}
}

pub struct HalfedgeIter<'a, M: Mesh> {
	id: HalfedgeId,
	mesh: &'a M,
}

impl <'a, M: Mesh> HalfedgeIter<'a, M> {
	pub fn new(id: HalfedgeId, mesh: &'a M) -> Self {
		Self { id, mesh }
	}

	fn len(&self) -> usize {
		self.mesh.n_halfedges()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_halfedges_capacity()
	}
}

impl <'a, M: Mesh> Deref for HalfedgeIter<'a, M> {
	type Target = HalfedgeId;

	fn deref(&self) -> &Self::Target {
		&self.id
	}
}

impl <'a, M: Mesh> Element for HalfedgeIter<'a, M> {
	type Item = HalfedgeId;
	type M = M;


	fn id(&self) -> HalfedgeId {
		self.id
	}

	fn mesh(&self) -> &M {
		self.mesh
	}

	fn valid(&self) -> bool {
		self.mesh.halfedge_is_valid(self.id)
	}

	fn next(&mut self) -> bool {
		self.id.0 += 1;
		self.id.0 < self.capacity()
	}
}

pub trait Halfedge: Element {
	fn edge(&self) -> EdgeIter<'_, Self::M>;
}

impl <'a, M: Mesh> Halfedge for HalfedgeIter<'a, M> {
	fn edge(&self) -> EdgeIter<'_, Self::M> {
		EdgeIter::new(self.mesh.he_edge(self.id), self.mesh)
	}
}

impl <'a, M: Mesh> Iterator for HalfedgeIter<'a, M> {
	type Item = HalfedgeId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}
