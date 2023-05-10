use bumpalo::{collections::Vec, vec, Bump};

use super::{
    abs_max, double_to_sign, dummy_asb_max, predicates, sign_reverse, GenericNum, ImplicitPoint3D,
    ImplicitPointLPI, ImplicitPointTPI, Orientation, Point3D,
};

/// Computes the orientation of the 3D points `pa`, `pb`, `pc` and `pd`.
/// If `pd` is above the plane defined by `pa`, `pb` and `pc`, then the
/// orientation is negative. If `pd` is below the plane, then the
/// orientation is positive. If `pd` is coplanar with `pa`, `pb` and `pc`,
/// then the orientation is zero.
pub fn orient2d_xy<'a, 'b: 'a>(
    pa: &'a Point3D<'b>,
    pb: &'a Point3D<'b>,
    pc: &'a Point3D<'b>,
    bump: &'b Bump,
) -> Orientation {
    match (pa, pb, pc) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            double_to_sign(predicates::orient2d(pa, pb, pc, bump))
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lee(pc, pa.data[0], pa.data[1], pb.data[0], pb.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tee(pc, pa.data[0], pa.data[1], pb.data[0], pb.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lee(pb, pc.data[0], pc.data[1], pa.data[0], pa.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => {
            orient2d_lle(pb, pc, pa.data[0], pa.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => {
            orient2d_lte(pb, pc, pa.data[0], pa.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tee(pb, pc.data[0], pc.data[1], pa.data[0], pa.data[1], bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => {
            sign_reverse(orient2d_lte(pc, pb, pa.data[0], pa.data[1], bump))
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => {
            orient2d_tte(pb, pc, pa.data[0], pa.data[1], bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_lee(pa, pb.data[0], pb.data[1], pc.data[0], pc.data[1], bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lle(pc, pa, pb.data[0], pb.data[1], bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            sign_reverse(orient2d_lte(pa, pc, pb.data[0], pb.data[1], bump))
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lle(pa, pb, pc.data[0], pc.data[1], bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte(pa, pb, pc.data[0], pc.data[1], bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt(pc, pa, pb),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt(pa, pb, pc),
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_tee(pa, pb.data[0], pb.data[1], pc.data[0], pc.data[1], bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lte(pc, pa, pb.data[0], pb.data[1], bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tte(pc, pa, pb.data[0], pb.data[1], bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            sign_reverse(orient2d_lte(pb, pa, pc.data[0], pc.data[1], bump))
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte(pa, pb, pc.data[0], pc.data[1], bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt(pc, pa, pb),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt(pa, pb, pc),
    }
}

fn orient2d_iee_impl<'a, 'b: 'a, const NEED_MAX: bool, T: 'b + GenericNum, F>(
    l1x: &'a T,
    l1y: &'a T,
    d1: &'a T,
    p2x: T,
    p2y: T,
    p3x: T,
    p3y: T,
    abs_max: F,
    bump: &'b Bump,
) -> (T, Option<T>)
where
    F: FnOnce(Vec<'b, T>) -> Option<T>,
{
    let t1x = &p2y - &p3y;
    let t1y = &p3x - &p2x;
    let e2 = l1x * &t1x;
    let e3 = l1y * &t1y;
    let e = e2 + e3;
    let pr1 = &p2x * &p3y;
    let pr2 = &p2y * &p3x;
    let pr = pr1 - pr2;
    let dpr = d1 * pr;
    let det = dpr + e;

    let max_var = if NEED_MAX {
        abs_max(vec![in bump; p2x, p2y, p3x, p3y, t1x, t1y])
    } else {
        None
    };
    (det, max_var)
}

fn orient2d_lee<'b>(
    p1: &ImplicitPointLPI,
    p2x: f64,
    p2y: f64,
    p3x: f64,
    p3y: f64,
    bump: &'b Bump,
) -> Orientation {
    if let Some(p1_static) = p1.static_filter() {
        let ret = orient2d_iee_impl::<'_, '_, true, _, _>(
            &p1_static.0.x,
            &p1_static.0.y,
            &p1_static.0.d,
            p2x,
            p2y,
            p3x,
            p3y,
            abs_max,
            bump,
        );
        let max_var = ret.1.unwrap().max(p1_static.1);
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= 4.752773695437811e-14;
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(p1_dynamic) = p1.dynamic_filter() {
        let (det, _) = orient2d_iee_impl::<'_, '_, false, _, _>(
            &p1_dynamic.x,
            &p1_dynamic.y,
            &p1_dynamic.d,
            p2x.into(),
            p2y.into(),
            p3x.into(),
            p3y.into(),
            dummy_asb_max,
            bump,
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(p1_exact) = p1.exact() {
        let (det, _) = orient2d_iee_impl::<'_, '_, false, _, _>(
            &p1_exact.x,
            &p1_exact.y,
            &p1_exact.d,
            vec![in bump; p2x].into(),
            vec![in bump; p2y].into(),
            vec![in bump; p3x].into(),
            vec![in bump; p3y].into(),
            dummy_asb_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient2d_tee<'b>(
    p1: &ImplicitPointTPI,
    p2x: f64,
    p2y: f64,
    p3x: f64,
    p3y: f64,
    bump: &'b Bump,
) -> Orientation {
    if let Some(p1_static) = p1.static_filter() {
        let ret = orient2d_iee_impl::<'_, '_, true, _, _>(
            &p1_static.0.x,
            &p1_static.0.y,
            &p1_static.0.d,
            p2x,
            p2y,
            p3x,
            p3y,
            abs_max,
            bump,
        );
        let max_var = ret.1.unwrap().max(p1_static.1);
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= 9.061883188277186e-13;
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(p1_dynamic) = p1.dynamic_filter() {
        let (det, _) = orient2d_iee_impl::<'_, '_, false, _, _>(
            &p1_dynamic.x,
            &p1_dynamic.y,
            &p1_dynamic.d,
            p2x.into(),
            p2y.into(),
            p3x.into(),
            p3y.into(),
            dummy_asb_max,
            bump,
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(p1_exact) = p1.exact() {
        let (det, _) = orient2d_iee_impl::<'_, '_, false, _, _>(
            &p1_exact.x,
            &p1_exact.y,
            &p1_exact.d,
            vec![in bump; p2x].into(),
            vec![in bump; p2y].into(),
            vec![in bump; p3x].into(),
            vec![in bump; p3y].into(),
            dummy_asb_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

fn orient2d_iie_impl<'a, 'b: 'a, const NEED_MAX: bool, T: 'b + GenericNum, F>(
    l1x: &'a T,
    l1y: &'a T,
    d1: &'a T,
    l2x: &'a T,
    l2y: &'a T,
    d2: &'a T,
    op3x: T,
    op3y: T,
    abs_max: F,
    bump: &'b Bump,
) -> (T, Option<T>)
where
    F: FnOnce(Vec<'b, T>) -> Option<T>,
{
    let a = d1 * l2x;
    let b = d2 * l1x;
    let c = d1 * &op3y;
    let e = d1 * l2y;
    let f = d2 * l1y;
    let g = d1 * &op3x;
    let ab = a - b;
    let cd = c - l1y;
    let ef = e - f;
    let gh = g - l1x;
    let abcd = ab * cd;
    let efgh = ef * gh;
    let det = abcd - efgh;
    let max_var = if NEED_MAX {
        abs_max(vec![in bump; op3x, op3y])
    } else {
        None
    };
    (det, max_var)
}

fn orient2d_iie<'b, IP1: ImplicitPoint3D<'b>, IP2: ImplicitPoint3D<'b>, F>(
    p1: &IP1,
    p2: &IP2,
    p3x: f64,
    p3y: f64,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation
where
    F: FnOnce(f64) -> f64,
{
    if let Some((p1_static, p2_static)) = p1.static_filter().zip(p2.static_filter()) {
        let ret = orient2d_iie_impl::<'_, '_, true, _, _>(
            &p1_static.0.x,
            &p1_static.0.y,
            &p1_static.0.d,
            &p2_static.0.x,
            &p2_static.0.y,
            &p2_static.0.d,
            p3x,
            p3y,
            abs_max,
            bump,
        );
        let max_var = ret.1.unwrap().max(p1_static.1).max(p2_static.1);
        let epsilon = static_filter_func(max_var);
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some((p1_dynamic, p2_dynamic)) = p1.dynamic_filter().zip(p2.dynamic_filter()) {
        let (det, _) = orient2d_iie_impl::<'_, '_, false, _, _>(
            &p1_dynamic.x,
            &p1_dynamic.y,
            &p1_dynamic.d,
            &p2_dynamic.x,
            &p2_dynamic.y,
            &p2_dynamic.d,
            p3x.into(),
            p3y.into(),
            dummy_asb_max,
            bump,
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some((p1_exact, p2_exact)) = p1.exact().zip(p2.exact()) {
        let (det, _) = orient2d_iie_impl::<'_, '_, false, _, _>(
            &p1_exact.x,
            &p1_exact.y,
            &p1_exact.d,
            &p2_exact.x,
            &p2_exact.y,
            &p2_exact.d,
            vec![in bump; p3x].into(),
            vec![in bump; p3y].into(),
            dummy_asb_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn orient2d_lle<'a, 'b: 'a>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointLPI<'b>,
    p3x: f64,
    p3y: f64,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie(
        p1,
        p2,
        p3x,
        p3y,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 1.699690735379461e-11
        },
        bump,
    )
}

#[inline(always)]
fn orient2d_lte<'a, 'b: 'a>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3x: f64,
    p3y: f64,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie(
        p1,
        p2,
        p3x,
        p3y,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 2.184958117212875e-10
        },
        bump,
    )
}

#[inline(always)]
fn orient2d_tte<'a, 'b: 'a>(
    p1: &'a ImplicitPointTPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3x: f64,
    p3y: f64,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie(
        p1,
        p2,
        p3x,
        p3y,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 3.307187945722514e-08
        },
        bump,
    )
}

fn orient2d_iii_impl<T: GenericNum>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    d2: &T,
    l3x: &T,
    l3y: &T,
    d3: &T,
) -> T {
    let a = d1 * l2x;
    let b = d2 * l1x;
    let c = d1 * l3y;
    let d = d3 * l1y;
    let e = d1 * l2y;
    let f = d2 * l1y;
    let g = d1 * l3x;
    let h = d3 * l1x;
    let ab = a - b;
    let cd = c - d;
    let ef = e - f;
    let gh = g - h;
    let abcd = ab * cd;
    let efgh = ef * gh;
    let det = abcd - efgh;
    det
}

fn orient2d_iii<
    'b,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    IP3: ImplicitPoint3D<'b>,
    F,
>(
    p1: &IP1,
    p2: &IP2,
    p3: &IP3,
    static_filter_func: F,
) -> Orientation
where
    F: FnOnce(f64) -> f64,
{
    if let Some(((p1_static, p2_static), p3_static)) = p1
        .static_filter()
        .zip(p2.static_filter())
        .zip(p3.static_filter())
    {
        let det = orient2d_iii_impl(
            &p1_static.0.x,
            &p1_static.0.y,
            &p1_static.0.d,
            &p2_static.0.x,
            &p2_static.0.y,
            &p2_static.0.d,
            &p3_static.0.x,
            &p3_static.0.y,
            &p3_static.0.d,
        );
        let max_var = p1_static.1.max(p2_static.1).max(p3_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(((p1_dynamic, p2_dynamic), p3_dynamic)) = p1
        .dynamic_filter()
        .zip(p2.dynamic_filter())
        .zip(p3.dynamic_filter())
    {
        let det = orient2d_iii_impl(
            &p1_dynamic.x,
            &p1_dynamic.y,
            &p1_dynamic.d,
            &p2_dynamic.x,
            &p2_dynamic.y,
            &p2_dynamic.d,
            &p3_dynamic.x,
            &p3_dynamic.y,
            &p3_dynamic.d,
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(((p1_exact, p2_exact), p3_exact)) = p1.exact().zip(p2.exact()).zip(p3.exact()) {
        let det = orient2d_iii_impl(
            &p1_exact.x,
            &p1_exact.y,
            &p1_exact.d,
            &p2_exact.x,
            &p2_exact.y,
            &p2_exact.d,
            &p3_exact.x,
            &p3_exact.y,
            &p3_exact.d,
        );
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn orient2d_lll<'b>(
    p1: &ImplicitPointLPI<'b>,
    p2: &ImplicitPointLPI<'b>,
    p3: &ImplicitPointLPI<'b>,
) -> Orientation {
    orient2d_iii(p1, p2, p3, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon * 1.75634284893534e-10
    })
}

#[inline(always)]
fn orient2d_llt<'b>(
    p1: &ImplicitPointLPI<'b>,
    p2: &ImplicitPointLPI<'b>,
    p3: &ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii(p1, p2, p3, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon * 2.144556754402072e-09
    })
}

#[inline(always)]
fn orient2d_ltt<'b>(
    p1: &ImplicitPointLPI<'b>,
    p2: &ImplicitPointTPI<'b>,
    p3: &ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii(p1, p2, p3, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon * 2.535681042914479e-08
    })
}

#[inline(always)]
fn orient2d_ttt<'b>(
    p1: &ImplicitPointTPI<'b>,
    p2: &ImplicitPointTPI<'b>,
    p3: &ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii(p1, p2, p3, |max_var| {
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
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon * 3.103174776697445e-06
    })
}
