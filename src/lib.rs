#![feature(const_fn_floating_point_arithmetic)]
#![feature(float_next_up_down)]
#![feature(cell_leak)]
#![feature(trait_alias)]
#![feature(iter_partition_in_place)]
#![feature(allocator_api)]
#![feature(portable_simd)]
#![feature(let_chains)]

use itertools::Itertools;

pub mod disjoint_set;
pub mod graphcut;
pub mod math;
pub mod mesh;
pub mod polygonlization;
pub mod predicates;
pub mod triangle;

const INVALID_IND: usize = usize::MAX;

#[inline(always)]
fn point(points: &[f64], tid: usize) -> &[f64] {
    let start = tid * 3;
    &points[start..(start + 3)]
}

#[inline(always)]
fn face_area_2d(points: &[f64]) -> f64 {
    points
        .chunks(2)
        .circular_tuple_windows()
        .map(|(pa, pb)| pa[0] * pb[1] - pa[1] * pb[0])
        .sum::<f64>()
}
