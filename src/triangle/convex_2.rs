use std::alloc::Allocator;

use itertools::Itertools;

use crate::{point_2, predicates::orient2d};

use super::alternate_axes;

pub fn convex_2<A: Allocator + Copy>(points: &[f64], alloc: A) -> Vec<usize, A> {
    let n_points = points.len() >> 1;
    let mut hull = Vec::with_capacity_in(n_points, alloc);
    hull.extend(0..n_points);
    hull.sort_unstable_by(|i: &usize, j: &usize| {
        let pa = point_2(points, *i);
        let pb = point_2(points, *j);
        pa.partial_cmp(pb).unwrap()
    });
    hull.dedup_by(|&mut i, &mut j| {
        let pa = point_2(points, i);
        let pb = point_2(points, j);
        pa == pb
    });

    alternate_axes(points, &mut hull, 0);
    let hull_len = div_conq_recurse(points, &mut hull, 0, &mut Vec::new_in(alloc), alloc);
    hull.resize(hull_len, 0);
    hull
}

fn div_conq_recurse<A: Allocator + Copy>(
    points: &[f64],
    hull: &mut [usize],
    axis: usize,
    merge: &mut Vec<usize, A>,
    alloc: A,
) -> usize {
    if hull.len() < 3 {
        return hull.len();
    } else if hull.len() == 3 {
        let pa = point_2(points, hull[0]);
        let pb = point_2(points, hull[1]);
        let pc = point_2(points, hull[2]);
        if orient2d(pa, pb, pc, alloc) < 0.0 {
            hull.swap(1, 2);
        }
        return 3;
    }
    let divider = hull.len() >> 1;

    let (left, right) = hull.split_at_mut(divider);
    let left_len = div_conq_recurse(points, left, 1 - axis, merge, alloc);
    let right_len = div_conq_recurse(points, right, 1 - axis, merge, alloc);

    merge_hulls(
        points,
        &left[..left_len],
        &right[..right_len],
        axis,
        merge,
        alloc,
    );
    for (tar, &src) in hull.iter_mut().zip(merge.iter()) {
        *tar = src;
    }

    return merge.len();
}

fn merge_hulls<A: Allocator + Copy>(
    points: &[f64],
    left: &[usize],
    right: &[usize],
    axis: usize,
    merge: &mut Vec<usize, A>,
    alloc: A,
) {
    let comp = |i: &&usize, j: &&usize| {
        let pa = point_2(points, **i);
        let pb = point_2(points, **j);
        if axis == 0 {
            pa.partial_cmp(pb).unwrap()
        } else {
            [pa[1], pa[0]].partial_cmp(&[pb[1], pb[0]]).unwrap()
        }
    };

    let left_max = left.iter().position_max_by(comp).unwrap();
    let right_min = right.iter().position_min_by(comp).unwrap();

    let [left_below, right_below] = get_tangent(points, left, left_max, right, right_min, alloc);
    let [right_top, left_top] = get_tangent(points, right, right_min, left, left_max, alloc);

    merge.clear();
    let mut idx = right_below;
    loop {
        merge.push(right[idx]);
        if idx == right_top {
            break;
        }
        idx = (idx + 1) % right.len();
    }
    idx = left_top;
    loop {
        merge.push(left[idx]);
        if idx == left_below {
            break;
        }
        idx = (idx + 1) % left.len();
    }
}

fn get_tangent<A: Allocator + Copy>(
    points: &[f64],
    left: &[usize],
    left_start: usize,
    right: &[usize],
    right_start: usize,
    alloc: A,
) -> [usize; 2] {
    let mut li = left_start;
    let mut ri = right_start;

    let next = |i: usize, len: usize| (i + 1) % len;
    let prev = |i: usize, len: usize| (i + len - 1) % len;

    loop {
        if left.len() > 1 {
            let lj = prev(li, left.len());
            let pa = point_2(points, left[li]);
            let pb = point_2(points, left[lj]);
            let pc = point_2(points, right[ri]);
            if orient2d(pa, pb, pc, alloc) > 0.0 {
                li = lj;
                continue;
            }
        }

        if right.len() > 1 {
            let rj = next(ri, right.len());
            let pa = point_2(points, left[li]);
            let pb = point_2(points, right[ri]);
            let pc = point_2(points, right[rj]);
            if orient2d(pa, pb, pc, alloc) < 0.0 {
                ri = rj;
                continue;
            }
        }
        break;
    }
    [li, ri]
}
