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

pub trait ElementId {
    fn new() -> Self;
    fn valid(&self) -> bool;
}

pub trait Element {
    type Item: From<usize>;
    type M: Mesh;
    fn id(&self) -> Self::Item;
    fn mesh(&self) -> &Self::M;
    fn valid(&self) -> bool;
    fn next(&mut self);
    fn is_end(&self) -> bool;
}

fn iter_next<E: Element>(ele: &mut E) -> Option<<E as Element>::Item> {
    if ele.is_end() {
        return None;
    }
    let id = ele.id();
    loop {
        ele.next();
        if ele.is_end() || ele.valid() {
            break;
        }
    }
    Some(id)
}
