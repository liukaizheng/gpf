mod vertex;
mod halfedge;
mod face;
mod boundary_loop;

pub use vertex::*;
pub use halfedge::*;
pub use face::*;
pub use boundary_loop::*;

use super::Mesh;

pub trait Element {
	type Item: From<usize>;
	type M : Mesh;
	fn id(&self) -> Self::Item;
	fn mesh(&self) -> &Self::M;
	fn len(&self) -> usize;
	fn capacity(&self) -> usize;
	fn valid(&self) -> bool;
	fn next(&mut self) -> bool;
}

fn iter_next<E: Element>(ele: &mut E) -> Option<<E as Element>::Item> {
	loop {
		if !ele.next() {
			return None;
		}
		if ele.valid() {
			return Some(ele.id());
		}
	}
}

// impl <'a, M: Mesh, E : Element<'a, M>> Iterator for E {
// 	type Item = Self::Item;

// 	fn next(&mut self) -> Option<Self::Item> {
// 		let ele_id = self.id + 1;
// 		while ele_id < self.capacity() && self.build(ele_id).valid() {
// 			ele_id += 1;
// 		}
// 		if ele_id >= self.capacity() {
// 			None
// 		} else {
// 			Some(Self::Item::from(ele_id))
// 		}
// 	}
// }