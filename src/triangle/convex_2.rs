use std::alloc::Allocator;

use crate::point_2;

pub fn convex_2<A: Allocator + Copy>(points: &[f64], alloc: A) -> Vec<f64, A> {
    let n_points = points.len() >> 1;
    let mut sorted_pt_inds = Vec::with_capacity_in(n_points, alloc);
    sorted_pt_inds.extend(0..n_points);
    let comp1 = |i: &usize, j: &usize| {
        let pa = point_2(points, *i);
        let pb = point_2(points, *j);
        pa.partial_cmp(pb).unwrap()
    };
    sorted_pt_inds.sort_unstable_by(comp1);
    sorted_pt_inds.dedup_by(|&mut i, &mut j| {
        let pa = point_2(points, i);
        let pb = point_2(points, j);
        pa == pb
    });
    let result = Vec::new_in(alloc);
    result
}
