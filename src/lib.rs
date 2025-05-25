#![feature(cell_leak)]
#![feature(trait_alias)]
#![feature(iter_partition_in_place)]
#![feature(allocator_api)]
#![feature(portable_simd)]
#![feature(let_chains)]
#![feature(iter_array_chunks)]
#![feature(array_chunks)]

use std::cell;

use itertools::Itertools;

pub mod boolean3d;
pub mod geometry;
pub mod graphcut;
pub mod math;
pub mod mesh;
pub mod polygonlization;
pub mod predicates;
pub mod triangle;
pub mod utils;

pub struct Tolerance {
    pub dist_tol: f64,
    pub cos_tol: f64,
}

impl Tolerance {
    pub fn new(dist_tol: f64, cos_tol: f64) -> Self {
        Tolerance { dist_tol, cos_tol }
    }
}
impl Default for Tolerance {
    fn default() -> Self {
        Tolerance {
            dist_tol: 1e-6,
            cos_tol: 1e-6,
        }
    }
}

const INVALID_IND: usize = usize::MAX;

#[inline(always)]
pub fn point<const N: usize>(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * N;
    &points[start..(start + N)]
}

#[inline(always)]
fn point_3(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

#[inline(always)]
pub fn point_2(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 2;
    &points[start..(start + 2)]
}

#[inline(always)]
pub fn oriented_index(idx: usize, reversed: bool) -> usize {
    (idx << 1) | if reversed { 1 } else { 0 }
}

#[inline(always)]
pub fn strip_orientation(idx: usize) -> usize {
    idx >> 1
}

#[inline(always)]
pub fn twin_index(idx: usize) -> usize {
    idx ^ 1
}

#[inline(always)]
pub fn is_positive(idx: usize) -> bool {
    idx & 1 == 0
}

#[inline(always)]
pub fn is_negative(idx: usize) -> bool {
    idx & 1 == 1
}

#[inline(always)]
pub fn face_area_2d(points: &[f64]) -> f64 {
    points
        .chunks(2)
        .circular_tuple_windows()
        .map(|(pa, pb)| pa[0] * pb[1] - pa[1] * pb[0])
        .sum::<f64>()
}
