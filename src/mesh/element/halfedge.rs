use std::ops::Index;

use crate::mesh::Mesh;
use super::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Halfedge(usize);

impl Index<Halfedge> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: Halfedge) -> &Self::Output {
		&self[index.0]
	}
}

pub struct HalfedgeBorrow<'a, M: Mesh> {
	id: usize,
	mesh: &'a M,
}

impl <'a, M: Mesh> HalfedgeBorrow<'a, M> {
	pub fn new(id: usize, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element<'a, M> for HalfedgeBorrow<'a, M> {
	type Item = Halfedge;

	fn build(&self, id: usize) -> Self {
		Self::new(id, self.mesh)
	}

	fn get(&self) -> Halfedge {
		Halfedge(self.id)
	}

	fn id(&self) -> usize {
		self.id
	}

	fn mesh(&self) -> &'a M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_halfedges()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_halfedges_capacity()
	}

	fn valid(&self) -> bool {
		self.id < self.capacity() && self.mesh.halfedge_is_valid(self.id)
	}
}
