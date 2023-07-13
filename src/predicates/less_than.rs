use std::alloc::Allocator;

use super::{
    double_to_sign, ExplicitPoint3D, GenericNum, ImplicitPoint3D, ImplicitPointLPI,
    ImplicitPointTPI, Orientation, Point3D,
};

pub fn less_than_on_x<A: Allocator + Copy>(pa: &Point3D, pb: &Point3D, bump: A) -> Orientation {
    match (pa, pb) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb)) => less_than_ee::<0>(pa, pb),
        (Point3D::Explicit(pa), Point3D::LPI(pb)) => less_than_on_x_le(pb, pa, bump),
        (Point3D::Explicit(pa), Point3D::TPI(pb)) => less_than_on_x_te(pb, pa, bump),
        (Point3D::LPI(pa), Point3D::Explicit(pb)) => less_than_on_x_le(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb)) => less_than_on_x_ll(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb)) => less_than_on_x_lt(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb)) => less_than_on_x_te(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb)) => less_than_on_x_lt(pb, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb)) => less_than_on_x_tt(pa, pb, bump),
    }
}

pub fn less_than_on_y<A: Allocator + Copy>(pa: &Point3D, pb: &Point3D, bump: A) -> Orientation {
    match (pa, pb) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb)) => less_than_ee::<1>(pa, pb),
        (Point3D::Explicit(pa), Point3D::LPI(pb)) => less_than_on_y_le(pb, pa, bump),
        (Point3D::Explicit(pa), Point3D::TPI(pb)) => less_than_on_y_te(pb, pa, bump),
        (Point3D::LPI(pa), Point3D::Explicit(pb)) => less_than_on_y_le(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb)) => less_than_on_y_ll(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb)) => less_than_on_y_lt(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb)) => less_than_on_y_te(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb)) => less_than_on_y_lt(pb, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb)) => less_than_on_y_tt(pa, pb, bump),
    }
}

pub fn less_than_on_z<A: Allocator + Copy>(pa: &Point3D, pb: &Point3D, bump: A) -> Orientation {
    match (pa, pb) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb)) => less_than_ee::<2>(pa, pb),
        (Point3D::Explicit(pa), Point3D::LPI(pb)) => less_than_on_z_le(pb, pa, bump),
        (Point3D::Explicit(pa), Point3D::TPI(pb)) => less_than_on_z_te(pb, pa, bump),
        (Point3D::LPI(pa), Point3D::Explicit(pb)) => less_than_on_z_le(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb)) => less_than_on_z_ll(pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb)) => less_than_on_z_lt(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb)) => less_than_on_z_te(pa, pb, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb)) => less_than_on_z_lt(pb, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb)) => less_than_on_z_tt(pa, pb, bump),
    }
}

#[inline(always)]
fn less_than_ee<const AXIS: usize>(pa: &ExplicitPoint3D, pb: &ExplicitPoint3D) -> Orientation {
    let a = pa[AXIS];
    let b = pb[AXIS];
    if a < b {
        return Orientation::Negative;
    } else if a > b {
        return Orientation::Positive;
    } else {
        return Orientation::Zero;
    }
}

#[inline(always)]
fn less_than_on_x_le<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_le::<0, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_y_le<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_le::<1, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_z_le<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_le::<2, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_x_te<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_te::<0, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_y_te<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_te::<1, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_z_te<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_te::<2, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_x_ll<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    less_than_ll::<0, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_y_ll<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    less_than_ll::<1, _>(pa, pb, bump)
}
#[inline(always)]
fn less_than_on_z_ll<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    less_than_ll::<2, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_x_lt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_lt::<0, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_y_lt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_lt::<1, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_z_lt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_lt::<2, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_x_tt<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_tt::<0, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_y_tt<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_tt::<1, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_on_z_tt<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_tt::<2, _>(pa, pb, bump)
}

#[inline(always)]
fn less_than_le<const AXIS: usize, A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_ie::<AXIS, _, _, _>(
        pa,
        pb,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon * 1.932297637868842e-14
        },
        bump,
    )
}

#[inline(always)]
fn less_than_te<const AXIS: usize, A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    less_than_ie::<AXIS, _, _, _>(
        pa,
        pb,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 3.980270973924514e-13
        },
        bump,
    )
}

#[inline(always)]
fn less_than_ll<const AXIS: usize, A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    less_than_ii::<AXIS, _, _, _, _>(
        pa,
        pb,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 2.922887626377607e-13
        },
        bump,
    )
}

#[inline(always)]
fn less_than_lt<const AXIS: usize, A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_ii::<AXIS, _, _, _, _>(
        pa,
        pb,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 4.321380059346694e-12
        },
        bump,
    )
}

#[inline(always)]
fn less_than_tt<const AXIS: usize, A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    less_than_ii::<AXIS, _, _, _, _>(
        pa,
        pb,
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
            epsilon * 5.504141586953918e-11
        },
        bump,
    )
}

#[inline(always)]
fn less_than_ie_impl<T: GenericNum>(l1x: &T, d1: &T, bx: T) -> T {
    let dbx = bx * d1;
    l1x - dbx
}

fn less_than_ie<
    const AXIS: usize,
    IP: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    pa: &IP,
    pb: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(pa_static) = pa.static_filter() {
        let (det, mut max_var) = if AXIS == 0 {
            (
                less_than_ie_impl(&pa_static.0.x, &pa_static.0.d, pb[0]),
                pb[0].abs(),
            )
        } else if AXIS == 1 {
            (
                less_than_ie_impl(&pa_static.0.y, &pa_static.0.d, pb[1]),
                pb[1].abs(),
            )
        } else {
            (
                less_than_ie_impl(&pa_static.0.z, &pa_static.0.d, pb[2]),
                pb[2].abs(),
            )
        };
        max_var = max_var.max(pa_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.dynamic_filter() {
        let det = if AXIS == 0 {
            less_than_ie_impl(&pa_dynamic.x, &pa_dynamic.d, pb[0].into())
        } else if AXIS == 1 {
            less_than_ie_impl(&pa_dynamic.y, &pa_dynamic.d, pb[1].into())
        } else {
            less_than_ie_impl(&pa_dynamic.z, &pa_dynamic.d, pb[2].into())
        };

        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else {
                return Orientation::Negative;
            }
        }
    }

    if let Some(pa_exact) = pa.exact(bump) {
        let det = if AXIS == 0 {
            less_than_ie_impl(
                &pa_exact.x,
                &pa_exact.d,
                [pb[0].into()].to_vec_in(bump).into(),
            )
        } else if AXIS == 1 {
            less_than_ie_impl(
                &pa_exact.y,
                &pa_exact.d,
                [pb[1].into()].to_vec_in(bump).into(),
            )
        } else {
            less_than_ie_impl(
                &pa_exact.z,
                &pa_exact.d,
                [pb[2].into()].to_vec_in(bump).into(),
            )
        };
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn less_than_ii_impl<T: GenericNum>(l1x: &T, d1: &T, l2x: &T, d2: &T) -> T {
    let k1 = d2 * l1x;
    let k2 = d1 * l2x;
    k1 - k2
}

fn less_than_ii<
    const AXIS: usize,
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    pa: &IP1,
    pb: &IP2,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(pa_static) = pa.static_filter() && let Some(pb_static) = pb.static_filter(){
        let det = if AXIS == 0 {
            less_than_ii_impl(&pa_static.0.x, &pa_static.0.d, &pb_static.0.x, &pb_static.0.d)
        } else if AXIS == 1 {
            less_than_ii_impl(&pa_static.0.y, &pa_static.0.d, &pb_static.0.y, &pb_static.0.d)
        } else {
            less_than_ii_impl(&pa_static.0.z, &pa_static.0.d, &pb_static.0.z, &pb_static.0.d)
        };
        let max_var = pa_static.1.max(pb_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.dynamic_filter() && let Some(pb_dynamic) = pb.dynamic_filter(){
        let det = if AXIS == 0 {
            less_than_ii_impl(&pa_dynamic.x, &pa_dynamic.d, &pb_dynamic.x, &pb_dynamic.d)
        } else if AXIS == 1 {
            less_than_ii_impl(&pa_dynamic.y, &pa_dynamic.d, &pb_dynamic.x, &pb_dynamic.d)
        } else {
            less_than_ii_impl(&pa_dynamic.z, &pa_dynamic.d, &pb_dynamic.x, &pb_dynamic.d)
        };

        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else {
                return Orientation::Negative;
            }
        }
    }

    if let Some(pa_exact) = pa.exact(bump) && let Some(pb_exact) = pb.exact(bump) {
        let det = if AXIS == 0 {
            less_than_ii_impl(
                &pa_exact.x,
                &pa_exact.d,
                &pb_exact.x,
                &pb_exact.d,
            )
        } else if AXIS == 1 {
            less_than_ii_impl(
                &pa_exact.y,
                &pa_exact.d,
                &pb_exact.y,
                &pb_exact.d,
            )
        } else {
            less_than_ii_impl(
                &pa_exact.z,
                &pa_exact.d,
                &pb_exact.z,
                &pb_exact.d,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}
