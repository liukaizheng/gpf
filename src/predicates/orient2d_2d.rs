use std::alloc::Allocator;

use crate::point;

use super::{
    abs_max, double_to_sign, dummy_abs_max,
    generic_point_2d::{ImplicitPointSSI, Point2D},
    GenericNum, Orientation,
};

pub fn orient2d<A: Allocator + Copy>(
    va: &Point2D,
    vb: &Point2D,
    vc: &Point2D,
    points: &[f64],
    alloc: A,
) -> Orientation {
    match (va, vb, vc) {
        (&Point2D::E(va), &Point2D::E(vb), &Point2D::E(vc)) => {
            let pa = point::<2>(points, va);
            let pb = point::<2>(points, vb);
            let pc = point::<2>(points, vc);
            double_to_sign(super::orient2d(pa, pb, pc, alloc))
        }
        (&Point2D::E(va), &Point2D::E(vb), Point2D::I(vc)) => {
            orient2d_iee(vc, va, vb, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(vb), &Point2D::E(vc)) => {
            orient2d_iee(vb, vc, va, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(vb), Point2D::I(vc)) => {
            orient2d_iie(vb, vc, va, points, alloc)
        }
        (Point2D::I(va), &Point2D::E(vb), &Point2D::E(vc)) => {
            orient2d_iee(va, vb, vc, points, alloc)
        }
        (Point2D::I(va), &Point2D::E(vb), Point2D::I(vc)) => {
            orient2d_iie(vc, va, vb, points, alloc)
        }
        (Point2D::I(va), Point2D::I(vb), &Point2D::E(vc)) => {
            orient2d_iie(va, vb, vc, points, alloc)
        }
        (Point2D::I(va), Point2D::I(vb), Point2D::I(vc)) => orient2d_iii(va, vb, vc, points, alloc),
    }
}

fn orient2d_iee_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    b1x: T,
    b1y: T,
    p2x: T,
    p2y: T,
    p3x: T,
    p3y: T,
    abs_max: F,
) -> (T, Option<T>) {
    let b1p3x = b1x - &p3x;
    let b1p3y = b1y - &p3y;
    let d1_b1p3x = d1 * &b1p3x;
    let d1_b1p3y = d1 * &b1p3y;
    let ix = d1_b1p3x + l1x;
    let iy = d1_b1p3y + l1y;
    let p2p3x = p2x - p3x;
    let p2p3y = p2y - p3y;
    let t0 = ix * &p2p3y;
    let t1 = iy * &p2p3x;
    let det = t0 - t1;

    if NEED_MAX {
        (det, abs_max(&[b1p3x, b1p3y, p2p3x, p2p3y]))
    } else {
        (det, None)
    }
}

fn orient2d_iee<A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    vb: usize,
    vc: usize,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb = point::<2>(points, vb);
    let pc = point::<2>(points, vc);
    if let Some(pa_static) = pa.ss_filter(points) {
        let ret = orient2d_iee_impl::<true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.d,
            pa_pa[0],
            pa_pa[1],
            pb[0],
            pb[1],
            pc[0],
            pc[1],
            abs_max,
        );
        let mut epsilon = ret.1.unwrap().max(pa_static.1);
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= 8.881784197001255e-15;
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points) {
        let ret = orient2d_iee_impl::<false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.d,
            pa_pa[0].into(),
            pa_pa[1].into(),
            pb[0].into(),
            pb[1].into(),
            pc[0].into(),
            pc[1].into(),
            dummy_abs_max,
        );

        if ret.0.positive() {
            return Orientation::Positive;
        } else if ret.0.negative() {
            return Orientation::Negative;
        }
    }

    if let Some(pa_exact) = pa.exact(points, alloc) {
        let (det, _) = orient2d_iee_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.d,
            [pa_pa[0]].to_vec_in(alloc).into(),
            [pa_pa[1]].to_vec_in(alloc).into(),
            [pb[0]].to_vec_in(alloc).into(),
            [pb[1]].to_vec_in(alloc).into(),
            [pc[0]].to_vec_in(alloc).into(),
            [pc[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient2d_iie_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    d2: &T,
    b1x: T,
    b1y: T,
    b2x: T,
    b2y: T,
    p3x: T,
    p3y: T,
    abs_max: F,
) -> (T, Option<T>) {
    let b1p3x = &b1x - &p3x;
    let b1p3y = &b1y - &p3y;
    let b2p3x = &b2x - &p3x;
    let b2p3y = &b2y - &p3y;
    let d1_b1p3x = d1 * &b1p3x;
    let d1_b1p3y = d1 * &b1p3y;
    let i1x = d1_b1p3x + l1x;
    let i1y = d1_b1p3y + l1y;
    let d2_b2p3x = d2 * &b2p3x;
    let d2_b2p3y = d2 * &b2p3y;
    let i2x = d2_b2p3x + l2x;
    let i2y = d2_b2p3y + l2y;
    let t0 = i1x * i2y;
    let t1 = i1y * i2x;
    let det = t0 - t1;

    if NEED_MAX {
        (det, abs_max(&[b1p3x, b1p3y, b2p3x, b2p3y]))
    } else {
        (det, None)
    }
}

fn orient2d_iie<A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    pb: &ImplicitPointSSI,
    vc: usize,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb_pa = point::<2>(points, pb.data[0]);
    let pc = point::<2>(points, vc);
    if let Some(pa_static) = pa.ss_filter(points)
        && let Some(pb_static) = pb.ss_filter(points)
    {
        let ret = orient2d_iie_impl::<true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.d,
            pa_pa[0],
            pa_pa[1],
            pb_pa[0],
            pb_pa[1],
            pc[0],
            pc[1],
            abs_max,
        );
        let max_var = ret.1.unwrap().max(pa_static.1).max(pb_static.1);
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= 5.684341886080809e-14;

        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points)
        && let Some(pb_dynamic) = pb.d_filter(points)
    {
        let ret = orient2d_iie_impl::<false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.d,
            pa_pa[0].into(),
            pa_pa[1].into(),
            pb_pa[0].into(),
            pb_pa[1].into(),
            pc[0].into(),
            pc[1].into(),
            dummy_abs_max,
        );

        if ret.0.positive() {
            return Orientation::Positive;
        } else if ret.0.negative() {
            return Orientation::Negative;
        }
    }

    if let Some(pa_exact) = pa.exact(points, alloc)
        && let Some(pb_exact) = pb.exact(points, alloc)
    {
        let (det, _) = orient2d_iie_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.d,
            [pa_pa[0]].to_vec_in(alloc).into(),
            [pa_pa[1]].to_vec_in(alloc).into(),
            [pb_pa[0]].to_vec_in(alloc).into(),
            [pb_pa[1]].to_vec_in(alloc).into(),
            [pc[0]].to_vec_in(alloc).into(),
            [pc[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient2d_iii_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    d2: &T,
    l3x: &T,
    l3y: &T,
    d3: &T,
    b1x: T,
    b1y: T,
    b2x: T,
    b2y: T,
    b3x: T,
    b3y: T,
    abs_max: F,
) -> (T, Option<T>) {
    let b1b3x = b1x - &b3x;
    let b1b3y = b1y - &b3y;
    let b2b3x = b2x - &b3x;
    let b2b3y = b2y - &b3y;
    let d1_b1b3x = d1 * &b1b3x;
    let d1_b1b3y = d1 * &b1b3y;
    let d1_b1b3_l1x = d1_b1b3x + l1x;
    let d1_b1b3_l1y = d1_b1b3y + l1y;
    let d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
    let d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
    let d2_b2b3x = d2 * &b2b3x;
    let d2_b2b3y = d2 * &b2b3y;
    let d2_b2b3_l2x = d2_b2b3x + l2x;
    let d2_b2b3_l2y = d2_b2b3y + l2y;
    let d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
    let d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
    let l3d1x = l3x * d1;
    let l3d1y = l3y * d1;
    let l3d2x = l3x * d2;
    let l3d2y = l3y * d2;
    let i1x = d3d1_b1b3_l1x - l3d1x;
    let i1y = d3d1_b1b3_l1y - l3d1y;
    let i2x = d3d2_b2b3_l2x - l3d2x;
    let i2y = d3d2_b2b3_l2y - l3d2y;
    let t0 = i1x * i2y;
    let t1 = i1y * i2x;
    let det = t0 - t1;

    if NEED_MAX {
        (det, abs_max(&[b1b3x, b1b3y, b2b3x, b2b3y]))
    } else {
        (det, None)
    }
}

fn orient2d_iii<'a, A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    pb: &ImplicitPointSSI,
    pc: &ImplicitPointSSI,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb_pa = point::<2>(points, pb.data[0]);
    let pc_pa = point::<2>(points, pc.data[0]);
    if let Some(pa_static) = pa.ss_filter(points)
        && let Some(pb_static) = pb.ss_filter(points)
        && let Some(pc_static) = pc.ss_filter(points)
    {
        let ret = orient2d_iii_impl::<true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.d,
            &pc_static.0.x,
            &pc_static.0.y,
            &pc_static.0.d,
            pa_pa[0],
            pa_pa[1],
            pb_pa[0],
            pb_pa[1],
            pc_pa[0],
            pc_pa[1],
            abs_max,
        );
        let max_var = ret.1.unwrap().max(pa_static.1).max(pb_static.1);
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= 5.684341886080809e-14;

        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points)
        && let Some(pb_dynamic) = pb.d_filter(points)
        && let Some(pc_dynamic) = pc.d_filter(points)
    {
        let ret = orient2d_iii_impl::<false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.d,
            &pc_dynamic.x,
            &pc_dynamic.y,
            &pc_dynamic.d,
            pa_pa[0].into(),
            pa_pa[1].into(),
            pb_pa[0].into(),
            pb_pa[1].into(),
            pc_pa[0].into(),
            pc_pa[1].into(),
            dummy_abs_max,
        );

        if ret.0.positive() {
            return Orientation::Positive;
        } else if ret.0.negative() {
            return Orientation::Negative;
        }
    }

    if let Some(pa_exact) = pa.exact(points, alloc)
        && let Some(pb_exact) = pb.exact(points, alloc)
        && let Some(pc_exact) = pc.exact(points, alloc)
    {
        let (det, _) = orient2d_iii_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.d,
            &pc_exact.x,
            &pc_exact.y,
            &pc_exact.d,
            [pa_pa[0]].to_vec_in(alloc).into(),
            [pa_pa[1]].to_vec_in(alloc).into(),
            [pb_pa[0]].to_vec_in(alloc).into(),
            [pb_pa[1]].to_vec_in(alloc).into(),
            [pc_pa[0]].to_vec_in(alloc).into(),
            [pc_pa[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}
