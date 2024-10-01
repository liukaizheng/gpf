use std::alloc::Allocator;

use crate::{point_2, predicates::orient2d};

use super::alternate_axes;

#[derive(Debug)]
struct LinkList<A: Allocator + Copy> {
    prev: Vec<usize, A>,
    next: Vec<usize, A>,
    real_indices: Vec<usize, A>,
}

impl<A: Allocator + Copy> LinkList<A> {
    fn new(n_points: usize, alloc: A) -> Self {
        let mut prev = Vec::with_capacity_in(n_points, alloc);
        let mut next = Vec::with_capacity_in(n_points, alloc);
        let mut real_indices = Vec::with_capacity_in(n_points, alloc);
        for i in 0..n_points {
            prev.push(i);
            next.push(i);
            real_indices.push(i);
        }
        LinkList {
            prev,
            next,
            real_indices,
        }
    }

    #[inline]
    fn prev(&self, i: usize) -> usize {
        self.prev[i]
    }

    #[inline]
    fn next(&self, i: usize) -> usize {
        self.next[i]
    }

    #[inline]
    fn connect(&mut self, i: usize, j: usize) {
        self.next[i] = j;
        self.prev[j] = i;
    }

    #[inline]
    fn duplicate(&mut self, i: usize) -> usize {
        let j = self.real_indices.len();
        let next = self.next[i];
        self.real_indices.push(i);
        self.prev.push(i);
        self.next.push(next);
        self.prev[next] = j;
        self.next[i] = j;
        j
    }

    #[inline]
    fn iter(&self, first: usize) -> LinkListIter<A> {
        LinkListIter {
            link: self,
            idx: first,
            first,
            finished: false,
        }
    }
}

struct LinkListIter<'a, A: Allocator + Copy> {
    link: &'a LinkList<A>,
    idx: usize,
    first: usize,
    finished: bool,
}

impl<'a, A: Allocator + Copy> Iterator for LinkListIter<'a, A> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let ret = self.idx;
        self.idx = self.link.next[self.idx];
        if self.idx == self.first {
            self.finished = true;
        }
        Some(ret)
    }
}

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
    let mut link = LinkList::new(n_points, alloc);
    let hull_start = div_conq_recurse(points, &hull, 0, &mut link, alloc);
    let mut result = Vec::new_in(alloc);
    result.extend(link.iter(hull_start).map(|idx| link.real_indices[idx]));
    result
}

fn div_conq_recurse<A: Allocator + Copy>(
    points: &[f64],
    hull: &[usize],
    axis: usize,
    link: &mut LinkList<A>,
    alloc: A,
) -> usize {
    match hull.len() {
        1 => {
            return hull[0];
        }
        2 => {
            link.connect(hull[0], hull[1]);
            link.connect(hull[1], hull[0]);
            return hull[0];
        }
        _ => {
            let divider = hull.len() >> 1;
            let (left, right) = hull.split_at(divider);
            let left_start = div_conq_recurse(points, left, 1 - axis, link, alloc);
            let right_start = div_conq_recurse(points, right, 1 - axis, link, alloc);

            return merge_hulls(points, left_start, right_start, axis, link, alloc);
        }
    }
}

fn merge_hulls<A: Allocator + Copy>(
    points: &[f64],
    left_start: usize,
    right_start: usize,
    axis: usize,
    link: &mut LinkList<A>,
    alloc: A,
) -> usize {
    let comp = |i: &usize, j: &usize| {
        let pa = point_2(points, link.real_indices[*i]);
        let pb = point_2(points, link.real_indices[*j]);
        if axis == 0 {
            pa.partial_cmp(pb).unwrap()
        } else {
            [pa[1], pa[0]].partial_cmp(&[pb[1], pb[0]]).unwrap()
        }
    };

    let left_max = link.iter(left_start).max_by(comp).unwrap();
    let right_min = link.iter(right_start).min_by(comp).unwrap();

    let [left_below, mut right_below] = get_tangent(points, left_max, right_min, link, alloc);
    let [right_top, mut left_top] = get_tangent(points, right_min, left_max, link, alloc);

    if left_below == left_top && right_below == right_top {
        if link.next(left_top) != left_top {
            left_top = link.duplicate(left_top);
        }
        if link.next(right_below) != right_below {
            right_below = link.duplicate(right_below);
        }
    }

    link.connect(left_below, right_below);
    link.connect(right_top, left_top);
    return left_below;
}

fn get_tangent<A: Allocator + Copy>(
    points: &[f64],
    left_start: usize,
    right_start: usize,
    link: &LinkList<A>,
    alloc: A,
) -> [usize; 2] {
    let mut left_next = left_start;
    let mut right_prev = right_start;

    loop {
        if link.prev(left_next) != left_next {
            let left_prev = link.prev(left_next);
            let pa = point_2(points, link.real_indices[left_next]);
            let pb = point_2(points, link.real_indices[left_prev]);
            let pc = point_2(points, link.real_indices[right_prev]);
            if orient2d(pa, pb, pc, alloc) > 0.0 {
                left_next = left_prev;
                continue;
            }
        }

        if link.next(right_prev) != right_prev {
            let right_next = link.next(right_prev);
            let pa = point_2(points, link.real_indices[left_next]);
            let pb = point_2(points, link.real_indices[right_prev]);
            let pc = point_2(points, link.real_indices[right_next]);
            if orient2d(pa, pb, pc, alloc) < 0.0 {
                right_prev = right_next;
                continue;
            }
        }
        break;
    }
    [left_next, right_prev]
}
