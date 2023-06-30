#![feature(const_fn_floating_point_arithmetic)]
#![feature(float_next_up_down)]
#![feature(cell_leak)]
#![feature(trait_alias)]
#![feature(iter_partition_in_place)]
#![feature(allocator_api)]
#![feature(portable_simd)]
#![feature(let_chains)]

pub mod math;
pub mod mesh;
pub mod polygonlization;
pub mod predicates;
pub mod triangle;
pub mod graphcut;

const INVALID_IND: usize = usize::MAX;
