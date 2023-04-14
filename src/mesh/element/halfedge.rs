use std::ops::Index;

use crate::mesh::Mesh;
use super::{Element, iter_next};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(usize);

impl From<usize> for HalfedgeId {
	fn from(id: usize) -> Self {
		Self(id)
	}
}

impl Index<HalfedgeId> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: HalfedgeId) -> &Self::Output {
		&self[index.0]
	}
}

pub struct Halfedge<'a, M: Mesh> {
	id: HalfedgeId,
	mesh: &'a M,
}

impl <'a, M: Mesh> Halfedge<'a, M> {
	pub fn new(id: HalfedgeId, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element for Halfedge<'a, M> {
	type Item = HalfedgeId;
	type M = M;


	fn id(&self) -> HalfedgeId {
		self.id
	}

	fn mesh(&self) -> &M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_halfedges()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_halfedges_capacity()
	}

	fn valid(&self) -> bool {
		self.mesh.halfedge_is_valid(self.id)
	}

	fn next(&mut self) -> bool {
		self.id.0 += 1;
		self.id.0 < self.capacity()
	}
}

impl <'a, M: Mesh> Iterator for Halfedge<'a, M> {
	type Item = HalfedgeId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}
