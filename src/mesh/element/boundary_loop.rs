use std::ops::Index;

use crate::mesh::Mesh;
use super::{Element, iter_next};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoopId(usize);

impl From<usize> for BoundaryLoopId {
	fn from(id: usize) -> Self {
		Self(id)
	}
}

impl Index<BoundaryLoopId> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: BoundaryLoopId) -> &Self::Output {
		&self[index.0]
	}
}

pub struct BoundaryLoop<'a, M: Mesh> {
	id: BoundaryLoopId,
	mesh: &'a M,
}

impl <'a, M: Mesh> BoundaryLoop<'a, M> {
	pub fn new(id: BoundaryLoopId, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element for BoundaryLoop<'a, M> {
	type Item = BoundaryLoopId;
	type M = M;


	fn id(&self) -> BoundaryLoopId {
		self.id
	}

	fn mesh(&self) -> &M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_boundary_loop()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_boundary_loop_capacity()
	}

	fn valid(&self) -> bool {
		self.mesh.boundary_loop_is_valid(self.id)
	}

	fn next(&mut self) -> bool {
		self.id.0 += 1;
		self.id.0 < self.capacity()
	}
}

impl <'a, M: Mesh> Iterator for BoundaryLoop<'a, M> {
	type Item = BoundaryLoopId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}