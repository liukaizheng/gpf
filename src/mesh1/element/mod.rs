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
