use std::ops::{Index, Deref};

use crate::mesh::Mesh;
use super::{Element, iter_next};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(usize);

impl From<usize> for VertexId {
	fn from(id: usize) -> Self {
		Self(id)
	}
}

impl <T> Index<VertexId> for Vec<T> {
	type Output = T;
	fn index(&self, index: VertexId) -> &Self::Output {
		&self[index.0]
	}
}

pub struct VertexIter<'a, M: Mesh> {
	id: VertexId,
	mesh: &'a M,
}

impl <'a, M: Mesh> VertexIter<'a, M> {
	pub fn new(id: VertexId, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Deref for VertexIter<'a, M> {
	type Target = VertexId;

	fn deref(&self) -> &Self::Target {
		&self.id
	}
}

impl <'a, M: Mesh> Element for VertexIter<'a, M> {
	type Item = VertexId;
	type M = M;

	fn id(&self) -> VertexId {
		self.id
	}
	fn mesh(&self) -> &M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_vertices()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_vertices_capacity()
	}

	fn valid(&self) -> bool {
	    self.mesh.vertex_is_valid(self.id)
	}

	fn next(&mut self) -> bool {
		self.id.0 += 1;
		self.id.0 < self.capacity()
	}
}

impl <'a, M: Mesh> Iterator for VertexIter<'a, M> {
	type Item = VertexId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}
