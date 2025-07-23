#[macro_use]
mod macros;

mod edge;
mod face;
mod halfedge;
mod vertex;

use crate::INVALID_IND;
use std::alloc::Allocator;

pub use edge::*;
pub use face::*;
pub use halfedge::*;
pub use vertex::*;

pub trait ElementId: From<usize> + Default {
    fn index(&self) -> usize;

    #[inline(always)]
    fn valid(&self) -> bool {
        self.index() != INVALID_IND
    }
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
