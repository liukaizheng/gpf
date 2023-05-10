#![feature(const_fn_floating_point_arithmetic)]
#![feature(float_next_up_down)]
#![feature(cell_leak)]
#![feature(trait_alias)]

pub mod mesh;
pub mod predicates;

const INVALID_IND: usize = usize::MAX;
