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
use bumpalo::collections::Vec;

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
pub fn ele_ranges<'b, E: ElementId>(
    start: usize,
    len: usize,
    bump: &'b bumpalo::Bump,
) -> Vec<'b, E> {
    Vec::from_iter_in((start..(start + len)).map(|idx| idx.into()), bump)
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
