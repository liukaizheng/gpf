mod boundary_loop;
mod edge;
mod face;
mod halfedge;
mod vertex;

use std::alloc::Allocator;

pub use boundary_loop::*;
pub use edge::*;
pub use face::*;
pub use halfedge::*;
pub use vertex::*;

use crate::INVALID_IND;

#[macro_use]
mod macros;

pub trait ElementIndex {
    fn index(&self) -> usize;
}
pub trait ElementId: From<usize> + ElementIndex {
    #[inline(always)]
    fn valid(&self) -> bool {
        self.index() != INVALID_IND
    }
}
pub trait Element {
    type Item;
    fn item(&self) -> Self::Item;
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

fn iter_next<E: Element>(ele: &mut E) -> Option<<E as Element>::Item> {
    if ele.is_end() {
        return None;
    }
    let item = ele.item();
    loop {
        ele.next();
        if ele.is_end() || ele.valid() {
            break;
        }
    }
    Some(item)
}
