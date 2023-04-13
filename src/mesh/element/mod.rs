mod vertex;
mod halfedge;
mod face;
mod boundary_loop;

pub use vertex::*;
pub use halfedge::*;
pub use face::*;
pub use boundary_loop::*;

use super::Mesh;

pub trait Element<'a, M: Mesh> {
	type Item: From<usize>;
	fn build(&self, id: usize) -> Self;
	fn get(&self) -> Self::Item;
	fn id(&self) -> usize;
	fn mesh(&self) -> &'a M;
	fn len(&self) -> usize;
	fn capacity(&self) -> usize;
	fn valid(&self) -> bool;
}

impl <'a, M: Mesh, E : Element<'a, M>> Iterator for E {
	type Item = Self::Item;

	fn next(&mut self) -> Option<Self::Item> {
		let ele_id = self.id + 1;
		while ele_id < self.capacity() && self.build(ele_id).valid() {
			ele_id += 1;
		}
		if ele_id >= self.capacity() {
			None
		} else {
			Some(Self::Item::from(ele_id))
		}
	}
}