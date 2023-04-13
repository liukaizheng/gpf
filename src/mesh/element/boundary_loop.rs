use std::ops::Index;

use crate::mesh::Mesh;
use super::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoop(usize);

impl Index<BoundaryLoop> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: BoundaryLoop) -> &Self::Output {
		&self[index.0]
	}
}

pub struct BoundaryLoopBorrow<'a, M: Mesh> {
	id: usize,
	mesh: &'a M,
}

impl <'a, M: Mesh> BoundaryLoopBorrow<'a, M> {
	pub fn new(id: usize, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element<'a, M> for BoundaryLoopBorrow<'a, M> {
	type Item = BoundaryLoop;

	fn build(&self, id: usize) -> Self {
		Self::new(id, self.mesh)
	}

	fn get(&self) -> BoundaryLoop {
		BoundaryLoop(self.id)
	}

	fn id(&self) -> usize {
		self.id
	}

	fn mesh(&self) -> &'a M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_boundary_loops()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_boundary_loops_capacity()
	}

	fn valid(&self) -> bool {
		self.id < self.capacity() && self.mesh.boundary_loop_is_valid(self.id)
	}
}