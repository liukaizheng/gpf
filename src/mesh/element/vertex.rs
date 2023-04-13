use std::ops::Index;

use crate::mesh::Mesh;
use super::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Vertex(usize);

impl Index<Vertex> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: Vertex) -> &Self::Output {
		&self[index.0]
	}
}

pub struct VertexBorrow<'a, M: Mesh> {
	id: usize,
	mesh: &'a M,
}

impl <'a, M: Mesh> VertexBorrow<'a, M> {
	pub fn new(id: usize, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element<'a, M> for VertexBorrow<'a, M> {
	type Item = Vertex;

	fn build(&self, id: usize) -> Self {
		Self::new(id, self.mesh)
	}

	fn get(&self) -> Vertex {
		Vertex(self.id)
	}

	fn id(&self) -> usize {
		self.id
	}
	fn mesh(&self) -> &'a M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_vertices()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_vertices_capacity()
	}

	fn valid(&self) -> bool {
		self.id < self.capacity() && self.mesh.vertex_is_valid(self.id)
	}
}
