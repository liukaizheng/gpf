use std::alloc::Allocator;

use crate::point;

use super::{
    abs_max, double_to_sign, dummy_abs_max,
    generic_point_2d::{ImplicitPointSSI, Point2D},
    predicates, GenericNum, Orientation,
};

pub fn incircle<T: AsRef<Point2D>, A: Allocator + Copy>(
    va: T,
    vb: T,
    vc: T,
    vd: T,
    points: &[f64],
    alloc: A,
) -> Orientation {
    match (va.as_ref(), vb.as_ref(), vc.as_ref(), vd.as_ref()) {
        (&Point2D::E(va), &Point2D::E(vb), &Point2D::E(vc), &Point2D::E(vd)) => {
            let pa = point::<2>(points, va);
            let pb = point::<2>(points, vb);
            let pc = point::<2>(points, vc);
            let pd = point::<2>(points, vd);
            double_to_sign(predicates::incircle(pa, pb, pc, pd, alloc))
        }
        (&Point2D::E(va), &Point2D::E(vb), &Point2D::E(vc), Point2D::I(pd)) => {
            incircle_ieee(pd, vb, va, vc, points, alloc)
        }
        (&Point2D::E(va), &Point2D::E(vb), Point2D::I(pc), &Point2D::E(vd)) => {
            incircle_ieee(pc, va, vb, vd, points, alloc)
        }
        (&Point2D::E(va), &Point2D::E(vb), Point2D::I(pc), Point2D::I(pd)) => {
            incircle_iiee(pc, pd, va, vb, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(pb), &Point2D::E(vc), &Point2D::E(vd)) => {
            incircle_ieee(pb, vc, va, vd, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(pb), &Point2D::E(vc), Point2D::I(pd)) => {
            incircle_iiee(pb, pd, vc, va, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(pb), Point2D::I(pc), &Point2D::E(vd)) => {
            incircle_iiee(pb, pc, va, vd, points, alloc)
        }
        (&Point2D::E(va), Point2D::I(pb), Point2D::I(pc), Point2D::I(pd)) => {
            incircle_iiie(pb, pd, pc, va, points, alloc)
        }
        (Point2D::I(pa), &Point2D::E(vb), &Point2D::E(vc), &Point2D::E(vd)) => {
            incircle_ieee(pa, vb, vc, vd, points, alloc)
        }
        (Point2D::I(pa), &Point2D::E(vb), &Point2D::E(vc), Point2D::I(pd)) => {
            incircle_iiee(pa, pd, vb, vc, points, alloc)
        }
        (Point2D::I(pa), &Point2D::E(vb), Point2D::I(pc), &Point2D::E(vd)) => {
            incircle_iiee(pa, pc, vd, vb, points, alloc)
        }
        (Point2D::I(pa), &Point2D::E(vb), Point2D::I(pc), Point2D::I(pd)) => {
            incircle_iiie(pa, pc, pd, vb, points, alloc)
        }
        (Point2D::I(pa), Point2D::I(pb), &Point2D::E(vc), &Point2D::E(vd)) => {
            incircle_iiee(pa, pb, vc, vd, points, alloc)
        }
        (Point2D::I(pa), Point2D::I(pb), &Point2D::E(vc), Point2D::I(pd)) => {
            incircle_iiie(pa, pd, pb, vc, points, alloc)
        }
        (Point2D::I(pa), Point2D::I(pb), Point2D::I(pc), &Point2D::E(vd)) => {
            incircle_iiie(pa, pb, pc, vd, points, alloc)
        }
        (Point2D::I(pa), Point2D::I(pb), Point2D::I(pc), Point2D::I(pd)) => {
            incircle_iiii(pa, pb, pc, pd, points, alloc)
        }
    }
}

fn incircle_ieee_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    b1x: T,
    b1y: T,
    pbx: T,
    pby: T,
    pcx: T,
    pcy: T,
    pdx: T,
    pdy: T,
    abs_max: F,
) -> (T, Option<T>) {
    let adx00 = b1x - &pdx;
    let adx0 = &adx00 * d1;
    let adx = adx0 + l1x;
    let ady00 = b1y - &pdy;
    let ady0 = &ady00 * d1;
    let ady = ady0 + l1y;
    let bdx = pbx - &pdx;
    let bdy = pby - &pdy;
    let cdx = pcx - pdx;
    let cdy = pcy - pdy;
    let abdeta = &adx * &bdy;
    let abdetb = &bdx * &ady;
    let abdet = abdeta - abdetb;
    let bcdeta = &bdx * &cdy;
    let bcdetb = &cdx * &bdy;
    let bcdet = bcdeta - bcdetb;
    let cadeta = &cdx * &ady;
    let cadetb = &adx * &cdy;
    let cadet = cadeta - cadetb;
    let alifta = &adx * &adx;
    let aliftb = &ady * &ady;
    let alift = alifta + aliftb;
    let blifta = &bdx * &bdx;
    let bliftb = &bdy * &bdy;
    let blift = blifta + bliftb;
    let clifta = &cdx * &cdx;
    let cliftb = &cdy * &cdy;
    let clift = clifta + cliftb;
    let la = alift * bcdet;
    let lb = blift * cadet;
    let lc = clift * abdet;
    let lbc = lb + lc;
    let lbcd = lbc * d1;
    let det = la + lbcd;
    if NEED_MAX {
        (det, abs_max(&[adx00, ady00, bdx, bdy, cdx, cdy]))
    } else {
        (det, None)
    }
}

fn incircle_ieee<A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    vb: usize,
    vc: usize,
    vd: usize,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb = point::<2>(points, vb);
    let pc = point::<2>(points, vc);
    let pd = point::<2>(points, vd);
    if let Some(pa_static) = pa.ss_filter(points) {
        let ret = incircle_ieee_impl::<true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.d,
            pa_pa[0],
            pa_pa[1],
            pb[0],
            pb[1],
            pc[0],
            pc[1],
            pd[0],
            pd[1],
            abs_max,
        );
        let mut epsilon = ret.1.unwrap().max(pa_static.1);

        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= 3.1263880373444464e-13;

        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points) {
        let ret = incircle_ieee_impl::<false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.d,
            pa_pa[0].into(),
            pa_pa[1].into(),
            pb[0].into(),
            pb[1].into(),
            pc[0].into(),
            pc[1].into(),
            pd[0].into(),
            pd[1].into(),
            dummy_abs_max,
        );

        if ret.0.positive() {
            return Orientation::Positive;
        } else if ret.0.negative() {
            return Orientation::Negative;
        }
    }

    if let Some(pa_exact) = pa.exact(points, alloc) {
        let (det, _) = incircle_ieee_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.d,
            [pa_pa[0]].to_vec_in(alloc).into(),
            [pa_pa[1]].to_vec_in(alloc).into(),
            [pb[0]].to_vec_in(alloc).into(),
            [pb[1]].to_vec_in(alloc).into(),
            [pc[0]].to_vec_in(alloc).into(),
            [pc[1]].to_vec_in(alloc).into(),
            [pd[0]].to_vec_in(alloc).into(),
            [pd[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn incircle_iiee_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
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
    pcx: T,
    pcy: T,
    pdx: T,
    pdy: T,
    abs_max: F,
) -> (T, Option<T>) {
    let beta_adx1 = b1x - &pdx;
    let beta_adx2 = &beta_adx1 * d1;
    let adx = beta_adx2 + l1x;
    let beta_ady1 = b1y - &pdy;
    let beta_ady2 = &beta_ady1 * d1;
    let ady = beta_ady2 + l1y;
    let beta_bdx1 = b2x - &pdx;
    let beta_bdx2 = &beta_bdx1 * d2;
    let bdx = beta_bdx2 + l2x;
    let beta_bdy1 = b2y - &pdy;
    let beta_bdy2 = &beta_bdy1 * d2;
    let bdy = beta_bdy2 + l2y;
    let cdx = pcx - pdx;
    let cdy = pcy - pdy;
    let abdeta = &adx * &bdy;
    let abdetb = &bdx * &ady;
    let abdet = abdeta - abdetb;
    let bcdeta = &bdx * &cdy;
    let bcdetb = &cdx * &bdy;
    let bcdet = bcdeta - bcdetb;
    let cadeta = &cdx * &ady;
    let cadetb = &adx * &cdy;
    let cadet = cadeta - cadetb;
    let alifta = &adx * &adx;
    let aliftb = &ady * &ady;
    let alift = alifta + aliftb;
    let blifta = &bdx * &bdx;
    let bliftb = &bdy * &bdy;
    let blift = blifta + bliftb;
    let clifta = &cdx * &cdx;
    let cliftb = &cdy * &cdy;
    let clift = clifta + cliftb;
    let la = alift * bcdet;
    let lb = blift * cadet;
    let lc = clift * abdet;
    let la1 = la * d2;
    let lb1 = lb * d1;
    let lc1 = lc * d1;
    let lc2 = lc1 * d2;
    let lab = la1 + lb1;
    let det = lab + lc2;

    if NEED_MAX {
        (
            det,
            abs_max(&[beta_adx1, beta_ady1, beta_bdx1, beta_bdy1, cdx, cdy]),
        )
    } else {
        (det, None)
    }
}

fn incircle_iiee<'a, A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    pb: &ImplicitPointSSI,
    vc: usize,
    vd: usize,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb_pa = point::<2>(points, pb.data[0]);
    let pc = point::<2>(points, vc);
    let pd = point::<2>(points, vd);
    if let Some(pa_static) = pa.ss_filter(points)
        && let Some(pb_static) = pb.ss_filter(points)
    {
        let ret = incircle_iiee_impl::<true, _, _>(
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
            pd[0],
            pd[1],
            abs_max,
        );
        let max_var = ret.1.unwrap().max(pa_static.1).max(pb_static.1);
        let mut epsilon = max_var;

        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= 4.746425474877485e-12;

        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points)
        && let Some(pb_dynamic) = pb.d_filter(points)
    {
        let ret = incircle_iiee_impl::<false, _, _>(
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
            pd[0].into(),
            pd[1].into(),
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
        let (det, _) = incircle_iiee_impl::<false, _, _>(
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
            [pd[0]].to_vec_in(alloc).into(),
            [pd[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}
fn incircle_iiie_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
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
    pdx: T,
    pdy: T,
    abs_max: F,
) -> (T, Option<T>) {
    let beta_adx1 = b1x - &pdx;
    let beta_adx2 = &beta_adx1 * d1;
    let adx = beta_adx2 + l1x;
    let beta_ady1 = b1y - &pdy;
    let beta_ady2 = &beta_ady1 * d1;
    let ady = beta_ady2 + l1y;
    let beta_bdx1 = b2x - &pdx;
    let beta_bdx2 = &beta_bdx1 * d2;
    let bdx = beta_bdx2 + l2x;
    let beta_bdy1 = b2y - &pdy;
    let beta_bdy2 = &beta_bdy1 * d2;
    let bdy = beta_bdy2 + l2y;
    let beta_cdx1 = b3x - &pdx;
    let beta_cdx2 = &beta_cdx1 * d3;
    let cdx = beta_cdx2 + l3x;
    let beta_cdy1 = b3y - pdy;
    let beta_cdy2 = &beta_cdy1 * d3;
    let cdy = beta_cdy2 + l3y;
    let abdeta = &adx * &bdy;
    let abdetb = &bdx * &ady;
    let abdet = abdeta - abdetb;
    let bcdeta = &bdx * &cdy;
    let bcdetb = &cdx * &bdy;
    let bcdet = bcdeta - bcdetb;
    let cadeta = &cdx * &ady;
    let cadetb = &adx * &cdy;
    let cadet = cadeta - cadetb;
    let alifta = &adx * &adx;
    let aliftb = &ady * &ady;
    let alift = alifta + aliftb;
    let blifta = &bdx * &bdx;
    let bliftb = &bdy * &bdy;
    let blift = blifta + bliftb;
    let clifta = &cdx * &cdx;
    let cliftb = &cdy * &cdy;
    let clift = clifta + cliftb;
    let la = alift * bcdet;
    let lb = blift * cadet;
    let lc = clift * abdet;
    let la1 = la * d2;
    let la2 = la1 * d3;
    let lb1 = lb * d1;
    let lb2 = lb1 * d3;
    let lc1 = lc * d1;
    let lc2 = lc1 * d2;
    let lab = la2 + lb2;
    let det = lab + lc2;
    if NEED_MAX {
        (
            det,
            abs_max(&[
                beta_adx1, beta_ady1, beta_bdx1, beta_bdy1, beta_cdx1, beta_cdy1,
            ]),
        )
    } else {
        (det, None)
    }
}

fn incircle_iiie<'a, A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    pb: &ImplicitPointSSI,
    pc: &ImplicitPointSSI,
    vd: usize,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb_pa = point::<2>(points, pb.data[0]);
    let pc_pa = point::<2>(points, pc.data[0]);
    let pd = point::<2>(points, vd);
    if let Some(pa_static) = pa.ss_filter(points)
        && let Some(pb_static) = pb.ss_filter(points)
        && let Some(pc_static) = pc.ss_filter(points)
    {
        let ret = incircle_iiie_impl::<true, _, _>(
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
            pd[0],
            pd[1],
            abs_max,
        );

        let max_var = ret.1.unwrap().max(pa_static.1).max(pb_static.1);
        let mut epsilon = max_var;

        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= 6.048139766790008e-11;
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
        let ret = incircle_iiie_impl::<false, _, _>(
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
            pd[0].into(),
            pd[1].into(),
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
        let (det, _) = incircle_iiie_impl::<false, _, _>(
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
            [pd[0]].to_vec_in(alloc).into(),
            [pd[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}
fn incircle_iiii_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    d2: &T,
    l3x: &T,
    l3y: &T,
    d3: &T,
    l4x: &T,
    l4y: &T,
    d4: &T,
    b1x: T,
    b1y: T,
    b2x: T,
    b2y: T,
    b3x: T,
    b3y: T,
    b4x: T,
    b4y: T,
    abs_max: F,
) -> (T, Option<T>) {
    let adx00 = l1x * d4;
    let adx01 = l4x * d1;
    let adx0 = adx00 - adx01;
    let adx10 = b1x - &b4x;
    let d1d4 = d1 * d4;
    let adx1 = &adx10 * &d1d4;
    let adx = adx0 + adx1;
    let ady00 = l1y * d4;
    let ady01 = l4y * d1;
    let ady0 = ady00 - ady01;
    let ady10 = b1y - &b4y;
    let ady1 = &ady10 * d1d4;
    let ady = ady0 + ady1;
    let bdx00 = l2x * d4;
    let bdx01 = l4x * d2;
    let bdx0 = bdx00 - bdx01;
    let bdx10 = b2x - &b4x;
    let d2d4 = d2 * d4;
    let bdx1 = &bdx10 * &d2d4;
    let bdx = bdx0 + bdx1;
    let bdy00 = l2y * d4;
    let bdy01 = l4y * d2;
    let bdy0 = bdy00 - bdy01;
    let bdy10 = b2y - &b4y;
    let bdy1 = &bdy10 * d2d4;
    let bdy = bdy0 + bdy1;
    let cdx00 = l3x * d4;
    let cdx01 = l4x * d3;
    let cdx0 = cdx00 - cdx01;
    let cdx10 = b3x - b4x;
    let d3d4 = d3 * d4;
    let cdx1 = &cdx10 * &d3d4;
    let cdx = cdx0 + cdx1;
    let cdy00 = l3y * d4;
    let cdy01 = l4y * d3;
    let cdy0 = cdy00 - cdy01;
    let cdy10 = b3y - b4y;
    let cdy1 = &cdy10 * d3d4;
    let cdy = cdy0 + cdy1;
    let abdeta = &adx * &bdy;
    let abdetb = &bdx * &ady;
    let abdet = abdeta - abdetb;
    let bcdeta = &bdx * &cdy;
    let bcdetb = &cdx * &bdy;
    let bcdet = bcdeta - bcdetb;
    let cadeta = &cdx * &ady;
    let cadetb = &adx * &cdy;
    let cadet = cadeta - cadetb;
    let alifta = &adx * &adx;
    let aliftb = &ady * &ady;
    let alift = alifta + aliftb;
    let blifta = &bdx * &bdx;
    let bliftb = &bdy * &bdy;
    let blift = blifta + bliftb;
    let clifta = &cdx * &cdx;
    let cliftb = &cdy * &cdy;
    let clift = clifta + cliftb;
    let la = alift * bcdet;
    let lb = blift * cadet;
    let lc = clift * abdet;
    let d2d3 = d2 * d3;
    let la1 = la * d2d3;
    let d1d3 = d1 * d3;
    let lb1 = lb * d1d3;
    let d1d2 = d1 * d2;
    let lc1 = lc * d1d2;
    let lab = la1 + lb1;
    let det = lab + lc1;
    if NEED_MAX {
        (det, abs_max(&[adx10, ady10, bdx10, bdy10, cdx10, cdy10]))
    } else {
        (det, None)
    }
}

fn incircle_iiii<'a, A: Allocator + Copy>(
    pa: &ImplicitPointSSI,
    pb: &ImplicitPointSSI,
    pc: &ImplicitPointSSI,
    pd: &ImplicitPointSSI,
    points: &[f64],
    alloc: A,
) -> Orientation {
    let pa_pa = point::<2>(points, pa.data[0]);
    let pb_pa = point::<2>(points, pb.data[0]);
    let pc_pa = point::<2>(points, pc.data[0]);
    let pd_pa = point::<2>(points, pd.data[0]);
    if let Some(pa_static) = pa.ss_filter(points)
        && let Some(pb_static) = pb.ss_filter(points)
        && let Some(pc_static) = pc.ss_filter(points)
        && let Some(pd_static) = pd.ss_filter(points)
    {
        let ret = incircle_iiii_impl::<true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.d,
            &pc_static.0.x,
            &pc_static.0.y,
            &pc_static.0.d,
            &pd_static.0.x,
            &pd_static.0.y,
            &pd_static.0.d,
            pa_pa[0],
            pa_pa[1],
            pb_pa[0],
            pb_pa[1],
            pc_pa[0],
            pc_pa[1],
            pd_pa[0],
            pd_pa[1],
            abs_max,
        );

        let max_var = ret.1.unwrap().max(pa_static.1).max(pb_static.1);
        let mut epsilon = max_var;

        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= 7.129983714548932e-09;

        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.d_filter(points)
        && let Some(pb_dynamic) = pb.d_filter(points)
        && let Some(pc_dynamic) = pc.d_filter(points)
        && let Some(pd_dynamic) = pd.d_filter(points)
    {
        let ret = incircle_iiii_impl::<false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.d,
            &pc_dynamic.x,
            &pc_dynamic.y,
            &pc_dynamic.d,
            &pd_dynamic.x,
            &pd_dynamic.y,
            &pd_dynamic.d,
            pa_pa[0].into(),
            pa_pa[1].into(),
            pb_pa[0].into(),
            pb_pa[1].into(),
            pc_pa[0].into(),
            pc_pa[1].into(),
            pd_pa[0].into(),
            pd_pa[1].into(),
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
        && let Some(pd_exact) = pd.exact(points, alloc)
    {
        let (det, _) = incircle_iiii_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.d,
            &pc_exact.x,
            &pc_exact.y,
            &pc_exact.d,
            &pd_exact.x,
            &pd_exact.y,
            &pd_exact.d,
            [pa_pa[0]].to_vec_in(alloc).into(),
            [pa_pa[1]].to_vec_in(alloc).into(),
            [pb_pa[0]].to_vec_in(alloc).into(),
            [pb_pa[1]].to_vec_in(alloc).into(),
            [pc_pa[0]].to_vec_in(alloc).into(),
            [pc_pa[1]].to_vec_in(alloc).into(),
            [pd_pa[0]].to_vec_in(alloc).into(),
            [pd_pa[1]].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}
