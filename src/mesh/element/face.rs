use std::ops::Index;

use crate::mesh::Mesh;
use super::{Element, iter_next};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(usize);

impl From<usize> for FaceId {
	fn from(id: usize) -> Self {
		Self(id)
	}
}

impl Index<FaceId> for Vec<usize> {
	type Output = usize;
	fn index(&self, index: FaceId) -> &Self::Output {
		&self[index.0]
	}
}

pub struct Face<'a, M: Mesh> {
	id: FaceId,
	mesh: &'a M,
}

impl <'a, M: Mesh> Face<'a, M> {
	pub fn new(id: FaceId, mesh: &'a M) -> Self {
		Self { id, mesh }
	}
}

impl <'a, M: Mesh> Element for Face<'a, M> {
	type Item = FaceId;
	type M = M;

	fn id(&self) -> FaceId {
		self.id
	}

	fn mesh(&self) -> &M {
		self.mesh
	}

	fn len(&self) -> usize {
		self.mesh.n_faces()
	}

	fn capacity(&self) -> usize {
		self.mesh.n_faces_capacity()
	}

	fn valid(&self) -> bool {
		self.mesh.face_is_valid(self.id)
	}

	fn next(&mut self) -> bool {
		self.id.0 += 1;
		self.id.0 < self.capacity()
	}
}

impl <'a, M: Mesh> Iterator for Face<'a, M> {
	type Item = FaceId;

	fn next(&mut self) -> Option<Self::Item> {
		iter_next(self)
	}
}