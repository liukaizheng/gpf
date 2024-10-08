use std::alloc::Allocator;

use crate::{
    point,
    predicates::{max_comp_in_tri_normal, orient2d, orient3d},
};

use super::triangulate;

pub enum Convex3Result<A: Allocator + Copy> {
    Dim0(usize),
    Dim1(Vec<usize, A>),
    Dim2(Vec<usize, A>),
    Dim3(Vec<usize, A>),
}

pub fn convex_3<A: Allocator + Copy>(points: &[f64], alloc: A) -> Convex3Result<A> {
    let n_points = points.len() / 3;
    let mut indices = Vec::with_capacity_in(n_points, alloc);
    indices.extend(0..n_points);

    indices.sort_unstable_by(|&i, &j| {
        let pa = point(points, i);
        let pb = point(points, j);
        pa.partial_cmp(pb).unwrap()
    });

    indices.dedup_by(|i, j| {
        let pa = point(points, *i);
        let pb = point(points, *j);
        pa[0] == pb[0] && pa[1] == pb[1] && pa[2] == pb[2]
    });

    if indices.len() == 1 {
        return Convex3Result::Dim0(indices[0]);
    }

    let colinear_len = hull_1(points, &indices, alloc);
    if colinear_len == indices.len() {
        return Convex3Result::Dim1(indices);
    }

    let coplanar_len = hull_2(points, &indices, colinear_len, alloc);

    let triangles = if coplanar_len == 3 {
        let mut triangles = Vec::with_capacity_in(3, alloc);
        triangles.extend(&indices[..3]);
        triangles
    } else {
        let axis = max_comp_in_tri_normal(
            point(points, indices[0]),
            point(points, indices[1]),
            point(points, indices[colinear_len]),
            alloc,
        );

        let mut points_2d = Vec::new_in(alloc);

        points_2d.extend(indices[..coplanar_len].into_iter().flat_map(|&idx| {
            let p = point(points, idx);
            [p[axis + 1], p[axis + 2]]
        }));
        triangulate(&points_2d, &[], alloc)
    };
    if colinear_len == indices.len() {
        return Convex3Result::Dim2(triangles);
    }
    return Convex3Result::Dim2(triangles);
}

fn hull_1<A: Allocator + Copy>(points: &[f64], indices: &[usize], alloc: A) -> usize {
    let pa = point(points, indices[0]);
    let pb = point(points, indices[1]);
    for i in 2..indices.len() {
        let pc = point(points, indices[i]);
        if orient2d(pa, pb, pc, alloc) != 0.0 {
            return i;
        }
    }
    return indices.len();
}

fn hull_2<A: Allocator + Copy>(points: &[f64], indices: &[usize], start: usize, alloc: A) -> usize {
    let pa = point(points, indices[0]);
    let pb = point(points, indices[1]);
    let pc = point(points, indices[start]);
    for i in (start + 1)..indices.len() {
        let pd = point(points, indices[i]);
        if orient3d(pa, pb, pc, pd, alloc) != 0.0 {
            return i;
        }
    }
    return indices.len();
}

fn hull_3<A: Allocator + Copy>(
    points: &[f64],
    indices: &[usize],
    start: usize,
    triangles: Vec<usize, A>,
    alloc: A,
) -> usize {
    0
}
