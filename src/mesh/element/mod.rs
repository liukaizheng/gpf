mod boundary_loop;
mod edge;
mod face;
mod halfedge;
mod vertex;

#[macro_use]
mod macros;

use std::alloc::Allocator;

pub use boundary_loop::*;
pub use edge::*;
pub use face::*;
pub use halfedge::*;
pub use vertex::*;

use super::Mesh;

pub trait ElementId: From<usize> {
    fn new() -> Self;
    fn valid(&self) -> bool;
}

pub trait Element {
    type Id: ElementId;
    type M: Mesh;
    fn id(&self) -> Self::Id;
    fn mesh(&self) -> &Self::M;
    fn valid(&self) -> bool;
    fn next(&mut self);
    fn is_end(&self) -> bool;
}

#[inline(always)]
pub fn ele_ranges<E: ElementId, A: Allocator + Copy>(
    start: usize,
    len: usize,
    bump: A,
) -> Vec<E, A> {
    let mut result = Vec::new_in(bump);
    result.extend((start..(start + len)).map(|idx| idx.into()));
    result
}

fn iter_next<E: Element>(ele: &mut E) -> Option<<E as Element>::Id> {
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
