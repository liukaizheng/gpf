use std::alloc::Allocator;

use super::{ExplicitPoint3D, GenericNum, ImplicitPoint3D, Orientation};

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
fn less_than_ie_impl<T: GenericNum>(l1x: &T, d1: &T, bx: T) -> T
where
{
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
