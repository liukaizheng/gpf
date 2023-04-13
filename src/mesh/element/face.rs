use std::ops::Index;

use crate::mesh::Mesh;
use super::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Face(usize);

impl Index<Face> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: Face) -> &Self::Output {
		&self[index.0]
	}
}

pub struct FaceBorrow<'a, M: Mesh> {
	id: usize,
	mesh: &'a M,
}

impl <'a, M: Mesh> FaceBorrow<'a, M> {
	pub fn new(id: usize, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element<'a, M> for FaceBorrow<'a, M> {
	type Item = Face;

	fn build(&self, id: usize) -> Self {
		Self::new(id, self.mesh)
	}

	fn get(&self) -> Face {
		Face(self.id)
	}

	fn id(&self) -> usize {
		self.id
	}

	fn mesh(&self) -> &'a M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_faces()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_faces_capacity()
	}

	fn valid(&self) -> bool {
		self.id < self.capacity() && self.mesh.face_is_valid(self.id)
	}
}
