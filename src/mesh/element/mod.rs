mod boundary_loop;
mod edge;
mod face;
mod halfedge;
mod vertex;

#[macro_use]
mod macros;

pub use boundary_loop::*;
pub use edge::*;
pub use face::*;
pub use halfedge::*;
pub use vertex::*;

use super::Mesh;

pub trait Element {
    type Item: From<usize>;
    type M: Mesh;
    fn id(&self) -> Self::Item;
    fn mesh(&self) -> &Self::M;
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
