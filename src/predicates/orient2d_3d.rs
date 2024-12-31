use std::alloc::Allocator;

use super::{
    abs_max, double_to_sign, dummy_abs_max, predicates, sign_reverse, ExplicitPoint3D, GenericNum,
    ImplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI, Orientation, Point3D,
};

#[inline(always)]
pub fn orient2d_by_axis<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    axis: usize,
    bump: A,
) -> Orientation {
    if axis == 0 {
        orient2d_yz(pa, pb, pc, bump)
    } else if axis == 1 {
        orient2d_zx(pa, pb, pc, bump)
    } else {
        orient2d_xy(pa, pb, pc, bump)
    }
}

/// Compute the orientation of the 3D points `pa`, `pb` and `pc`.
/// If `pc` is above the line defined by `pa` and `pb`, then the
/// orientation is negative. If `pc` is below the line, then the
/// orientation is positive. If `pc` is collinear with `pa` and `pb`,
/// then the orientation is zero.
pub fn orient2d_xy<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    bump: A,
) -> Orientation {
    match (pa, pb, pc) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            double_to_sign(predicates::orient2d(pa, pb, pc, bump))
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lee_xy(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tee_xy(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_xy(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => {
            orient2d_lle_xy(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => {
            orient2d_lte_xy(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_xy(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => {
            sign_reverse(orient2d_lte_xy(pc, pb, pa, bump))
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => {
            orient2d_tte_xy(pb, pc, pa, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_xy(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lle_xy(pc, pa, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            sign_reverse(orient2d_lte_xy(pa, pc, pb, bump))
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lle_xy(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_xy(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_xy(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_xy(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_xy(pc, pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_xy(pa, pb, pc, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_xy(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lte_xy(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tte_xy(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            sign_reverse(orient2d_lte_xy(pb, pa, pc, bump))
        }

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_xy(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_xy(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_xy(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_xy(pc, pa, pb, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_xy(pa, pb, pc, bump),
    }
}

pub fn orient2d_zx<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    bump: A,
) -> Orientation {
    match (pa, pb, pc) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            double_to_sign(predicates::orient2d(
                &[pa.data[2], pa.data[0]],
                &[pb.data[2], pb.data[0]],
                &[pc.data[2], pc.data[0]],
                bump,
            ))
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lee_zx(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tee_zx(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_zx(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => {
            orient2d_lle_zx(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => {
            orient2d_lte_zx(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_zx(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => {
            sign_reverse(orient2d_lte_zx(pc, pb, pa, bump))
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => {
            orient2d_tte_zx(pb, pc, pa, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_zx(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lle_zx(pc, pa, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            sign_reverse(orient2d_lte_zx(pa, pc, pb, bump))
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lle_zx(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_zx(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_zx(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_zx(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_zx(pc, pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_zx(pa, pb, pc, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_zx(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lte_zx(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tte_zx(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            sign_reverse(orient2d_lte_zx(pb, pa, pc, bump))
        }

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_zx(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_zx(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_zx(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_zx(pc, pa, pb, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_zx(pa, pb, pc, bump),
    }
}

pub fn orient2d_yz<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    bump: A,
) -> Orientation {
    match (pa, pb, pc) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            double_to_sign(predicates::orient2d(&pa[1..], &pb[1..], &pc[1..], bump))
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lee_yz(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tee_yz(pc, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_yz(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => {
            orient2d_lle_yz(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => {
            orient2d_lte_yz(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_yz(pb, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => {
            sign_reverse(orient2d_lte_yz(pc, pb, pa, bump))
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => {
            orient2d_tte_yz(pb, pc, pa, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_lee_yz(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lle_yz(pc, pa, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            sign_reverse(orient2d_lte_yz(pa, pc, pb, bump))
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lle_yz(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_yz(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_yz(pa, pb, pc, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_yz(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_yz(pc, pa, pb, bump),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_yz(pa, pb, pc, bump),
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            orient2d_tee_yz(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc)) => {
            orient2d_lte_yz(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc)) => {
            orient2d_tte_yz(pc, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc)) => {
            sign_reverse(orient2d_lte_yz(pb, pa, pc, bump))
        }

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_yz(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_yz(pb, pc, pa, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_yz(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_yz(pc, pa, pb, bump),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_yz(pa, pb, pc, bump),
    }
}

/// Implicit-Explicit-Explicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lee_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lee::<2, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_xy<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tee::<2, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lle_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lle::<2, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lte::<2, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_xy<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tte::<2, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lll_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    allocator: A,
) -> Orientation {
    orient2d_lll::<2, _>(pa, pb, pc, allocator)
}
#[inline(always)]
fn orient2d_llt_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_llt::<2, _>(pa, pb, pc, allocator)
}
#[inline(always)]
fn orient2d_ltt_xy<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_ltt::<2, _>(pa, pb, pc, allocator)
}
#[inline(always)]
fn orient2d_ttt_xy<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_ttt::<2, _>(pa, pb, pc, allocator)
}

/// Implicit-Explicit-Explicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lee_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lee::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_zx<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tee::<1, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lle_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lle::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lte::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_zx<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tte::<1, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lll_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    orient2d_lll::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_llt_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_llt::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_ltt_zx<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_ltt::<1, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_ttt_zx<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_ttt::<1, _>(pa, pb, pc, bump)
}

/// Implicit-Explicit-Explicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lee_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lee::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_yz<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tee::<0, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lle_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lle::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_lte::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_yz<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_tte::<0, _>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lll_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    orient2d_lll::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_llt_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_llt::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_ltt_yz<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_ltt::<0, _>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_ttt_yz<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient2d_ttt::<0, _>(pa, pb, pc, bump)
}

fn orient2d_iee_impl<const NEED_MAX: bool, T: GenericNum, F>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    p2x: T,
    p2y: T,
    p3x: T,
    p3y: T,
    abs_max: F,
) -> (T, Option<T>)
where
    F: FnOnce(&[T]) -> Option<T>,
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
        abs_max(&[p2x, p2y, p3x, p3y, t1x, t1y])
    } else {
        None
    };
    (det, max_var)
}
fn orient2d_iee<
    const AXIS: u32,
    IP: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    p1: &IP,
    p2: &ExplicitPoint3D,
    p3: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(p1_static) = p1.static_filter() {
        let ret = if AXIS == 2 {
            orient2d_iee_impl::<true, _, _>(
                &p1_static.0.x,
                &p1_static.0.y,
                &p1_static.0.d,
                p2.data[0],
                p2.data[1],
                p3.data[0],
                p3.data[1],
                abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<true, _, _>(
                &p1_static.0.z,
                &p1_static.0.x,
                &p1_static.0.d,
                p2.data[2],
                p2.data[0],
                p3.data[2],
                p3.data[0],
                abs_max,
            )
        } else {
            orient2d_iee_impl::<true, _, _>(
                &p1_static.0.y,
                &p1_static.0.z,
                &p1_static.0.d,
                p2.data[1],
                p2.data[2],
                p3.data[1],
                p3.data[2],
                abs_max,
            )
        };
        let max_var = ret.1.unwrap().max(p1_static.1);
        let epsilon = static_filter_func(max_var);
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(p1_dynamic) = p1.dynamic_filter() {
        let (det, _) = if AXIS == 2 {
            orient2d_iee_impl::<false, _, _>(
                &p1_dynamic.x,
                &p1_dynamic.y,
                &p1_dynamic.d,
                p2.data[0].into(),
                p2.data[1].into(),
                p3.data[0].into(),
                p3.data[1].into(),
                dummy_abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<false, _, _>(
                &p1_dynamic.z,
                &p1_dynamic.x,
                &p1_dynamic.d,
                p2.data[2].into(),
                p2.data[0].into(),
                p3.data[2].into(),
                p3.data[0].into(),
                dummy_abs_max,
            )
        } else {
            orient2d_iee_impl::<false, _, _>(
                &p1_dynamic.y,
                &p1_dynamic.z,
                &p1_dynamic.d,
                p2.data[1].into(),
                p2.data[2].into(),
                p3.data[1].into(),
                p3.data[2].into(),
                dummy_abs_max,
            )
        };
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(p1_exact) = p1.exact(bump) {
        let (det, _) = if AXIS == 2 {
            orient2d_iee_impl::<false, _, _>(
                &p1_exact.x,
                &p1_exact.y,
                &p1_exact.d,
                [p2.data[0]].to_vec_in(bump).into(),
                [p2.data[1]].to_vec_in(bump).into(),
                [p3.data[0]].to_vec_in(bump).into(),
                [p3.data[1]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<false, _, _>(
                &p1_exact.z,
                &p1_exact.x,
                &p1_exact.d,
                [p2.data[2]].to_vec_in(bump).into(),
                [p2.data[0]].to_vec_in(bump).into(),
                [p3.data[2]].to_vec_in(bump).into(),
                [p3.data[0]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        } else {
            orient2d_iee_impl::<false, _, _>(
                &p1_exact.y,
                &p1_exact.z,
                &p1_exact.d,
                [p2.data[1]].to_vec_in(bump).into(),
                [p2.data[2]].to_vec_in(bump).into(),
                [p3.data[1]].to_vec_in(bump).into(),
                [p3.data[2]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient2d_lee<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ExplicitPoint3D,
    p3: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_iee::<AXIS, _, _, _>(
        p1,
        p2,
        p3,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon * 4.752773695437811e-14
        },
        bump,
    )
}

fn orient2d_tee<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointTPI,
    p2: &ExplicitPoint3D,
    p3: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_iee::<AXIS, _, _, _>(
        p1,
        p2,
        p3,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon * 9.061883188277186e-13
        },
        bump,
    )
}

fn orient2d_iie_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    d2: &T,
    op3x: T,
    op3y: T,
    abs_max: F,
) -> (T, Option<T>) {
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
        abs_max(&[op3x, op3y])
    } else {
        None
    };
    (det, max_var)
}

fn orient2d_iie<
    const AXIS: u32,
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    A: Allocator + Copy,
    F: FnOnce(f64) -> f64,
>(
    p1: &IP1,
    p2: &IP2,
    p3: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some((p1_static, p2_static)) = p1.static_filter().zip(p2.static_filter()) {
        let ret = if AXIS == 2 {
            orient2d_iie_impl::<true, _, _>(
                &p1_static.0.x,
                &p1_static.0.y,
                &p1_static.0.d,
                &p2_static.0.x,
                &p2_static.0.y,
                &p2_static.0.d,
                p3.data[0],
                p3.data[1],
                abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<true, _, _>(
                &p1_static.0.z,
                &p1_static.0.x,
                &p1_static.0.d,
                &p2_static.0.z,
                &p2_static.0.x,
                &p2_static.0.d,
                p3.data[2],
                p3.data[0],
                abs_max,
            )
        } else {
            orient2d_iie_impl::<true, _, _>(
                &p1_static.0.y,
                &p1_static.0.z,
                &p1_static.0.d,
                &p2_static.0.y,
                &p2_static.0.z,
                &p2_static.0.d,
                p3.data[1],
                p3.data[2],
                abs_max,
            )
        };
        let max_var = ret.1.unwrap().max(p1_static.1).max(p2_static.1);
        let epsilon = static_filter_func(max_var);
        if ret.0 > epsilon {
            return Orientation::Positive;
        } else if ret.0 < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some((p1_dynamic, p2_dynamic)) = p1.dynamic_filter().zip(p2.dynamic_filter()) {
        let (det, _) = if AXIS == 2 {
            orient2d_iie_impl::<false, _, _>(
                &p1_dynamic.x,
                &p1_dynamic.y,
                &p1_dynamic.d,
                &p2_dynamic.x,
                &p2_dynamic.y,
                &p2_dynamic.d,
                p3.data[0].into(),
                p3.data[1].into(),
                dummy_abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<false, _, _>(
                &p1_dynamic.z,
                &p1_dynamic.x,
                &p1_dynamic.d,
                &p2_dynamic.z,
                &p2_dynamic.x,
                &p2_dynamic.d,
                p3.data[2].into(),
                p3.data[0].into(),
                dummy_abs_max,
            )
        } else {
            orient2d_iie_impl::<false, _, _>(
                &p1_dynamic.y,
                &p1_dynamic.z,
                &p1_dynamic.d,
                &p2_dynamic.y,
                &p2_dynamic.z,
                &p2_dynamic.d,
                p3.data[1].into(),
                p3.data[2].into(),
                dummy_abs_max,
            )
        };
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some((p1_exact, p2_exact)) = p1.exact(bump).zip(p2.exact(bump)) {
        let (det, _) = if AXIS == 2 {
            orient2d_iie_impl::<false, _, _>(
                &p1_exact.x,
                &p1_exact.y,
                &p1_exact.d,
                &p2_exact.x,
                &p2_exact.y,
                &p2_exact.d,
                [p3.data[0]].to_vec_in(bump).into(),
                [p3.data[1]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<false, _, _>(
                &p1_exact.z,
                &p1_exact.x,
                &p1_exact.d,
                &p2_exact.z,
                &p2_exact.x,
                &p2_exact.d,
                [p3.data[2]].to_vec_in(bump).into(),
                [p3.data[0]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        } else {
            orient2d_iie_impl::<false, _, _>(
                &p1_exact.y,
                &p1_exact.z,
                &p1_exact.d,
                &p2_exact.y,
                &p2_exact.z,
                &p2_exact.d,
                [p3.data[1]].to_vec_in(bump).into(),
                [p3.data[2]].to_vec_in(bump).into(),
                dummy_abs_max,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn orient2d_lle<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ImplicitPointLPI,
    p3: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_iie::<AXIS, _, _, _, _>(
        p1,
        p2,
        p3,
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
fn orient2d_lte<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ImplicitPointTPI,
    p3: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_iie::<AXIS, _, _, _, _>(
        p1,
        p2,
        p3,
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
fn orient2d_tte<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointTPI,
    p2: &ImplicitPointTPI,
    p3: &ExplicitPoint3D,
    bump: A,
) -> Orientation {
    orient2d_iie::<AXIS, _, _, _, _>(
        p1,
        p2,
        p3,
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
    const AXIS: u32,
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    IP3: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    p1: &IP1,
    p2: &IP2,
    p3: &IP3,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(((p1_static, p2_static), p3_static)) = p1
        .static_filter()
        .zip(p2.static_filter())
        .zip(p3.static_filter())
    {
        let det = if AXIS == 2 {
            orient2d_iii_impl(
                &p1_static.0.x,
                &p1_static.0.y,
                &p1_static.0.d,
                &p2_static.0.x,
                &p2_static.0.y,
                &p2_static.0.d,
                &p3_static.0.x,
                &p3_static.0.y,
                &p3_static.0.d,
            )
        } else if AXIS == 1 {
            orient2d_iii_impl(
                &p1_static.0.z,
                &p1_static.0.x,
                &p1_static.0.d,
                &p2_static.0.z,
                &p2_static.0.x,
                &p2_static.0.d,
                &p3_static.0.z,
                &p3_static.0.x,
                &p3_static.0.d,
            )
        } else {
            orient2d_iii_impl(
                &p1_static.0.y,
                &p1_static.0.z,
                &p1_static.0.d,
                &p2_static.0.y,
                &p2_static.0.z,
                &p2_static.0.d,
                &p3_static.0.y,
                &p3_static.0.z,
                &p3_static.0.d,
            )
        };
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
        let det = if AXIS == 2 {
            orient2d_iii_impl(
                &p1_dynamic.x,
                &p1_dynamic.y,
                &p1_dynamic.d,
                &p2_dynamic.x,
                &p2_dynamic.y,
                &p2_dynamic.d,
                &p3_dynamic.x,
                &p3_dynamic.y,
                &p3_dynamic.d,
            )
        } else if AXIS == 1 {
            orient2d_iii_impl(
                &p1_dynamic.z,
                &p1_dynamic.x,
                &p1_dynamic.d,
                &p2_dynamic.z,
                &p2_dynamic.x,
                &p2_dynamic.d,
                &p3_dynamic.z,
                &p3_dynamic.x,
                &p3_dynamic.d,
            )
        } else {
            orient2d_iii_impl(
                &p1_dynamic.y,
                &p1_dynamic.z,
                &p1_dynamic.d,
                &p2_dynamic.y,
                &p2_dynamic.z,
                &p2_dynamic.d,
                &p3_dynamic.y,
                &p3_dynamic.z,
                &p3_dynamic.d,
            )
        };
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(((p1_exact, p2_exact), p3_exact)) =
        p1.exact(bump).zip(p2.exact(bump)).zip(p3.exact(bump))
    {
        let det = if AXIS == 2 {
            orient2d_iii_impl(
                &p1_exact.x,
                &p1_exact.y,
                &p1_exact.d,
                &p2_exact.x,
                &p2_exact.y,
                &p2_exact.d,
                &p3_exact.x,
                &p3_exact.y,
                &p3_exact.d,
            )
        } else if AXIS == 1 {
            orient2d_iii_impl(
                &p1_exact.z,
                &p1_exact.x,
                &p1_exact.d,
                &p2_exact.z,
                &p2_exact.x,
                &p2_exact.d,
                &p3_exact.z,
                &p3_exact.x,
                &p3_exact.d,
            )
        } else {
            orient2d_iii_impl(
                &p1_exact.y,
                &p1_exact.z,
                &p1_exact.d,
                &p2_exact.y,
                &p2_exact.z,
                &p2_exact.d,
                &p3_exact.y,
                &p3_exact.z,
                &p3_exact.d,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn orient2d_lll<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ImplicitPointLPI,
    p3: &ImplicitPointLPI,
    allocator: A,
) -> Orientation {
    orient2d_iii::<AXIS, _, _, _, _, _>(
        p1,
        p2,
        p3,
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
            epsilon * 1.75634284893534e-10
        },
        allocator,
    )
}

#[inline(always)]
fn orient2d_llt<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ImplicitPointLPI,
    p3: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_iii::<AXIS, _, _, _, _, _>(
        p1,
        p2,
        p3,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon * 2.144556754402072e-09
        },
        allocator,
    )
}

#[inline(always)]
fn orient2d_ltt<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointLPI,
    p2: &ImplicitPointTPI,
    p3: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_iii::<AXIS, _, _, _, _, _>(
        p1,
        p2,
        p3,
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
            epsilon * 2.535681042914479e-08
        },
        allocator,
    )
}

#[inline(always)]
fn orient2d_ttt<const AXIS: u32, A: Allocator + Copy>(
    p1: &ImplicitPointTPI,
    p2: &ImplicitPointTPI,
    p3: &ImplicitPointTPI,
    allocator: A,
) -> Orientation {
    orient2d_iii::<AXIS, _, _, _, _, _>(
        p1,
        p2,
        p3,
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
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 3.103174776697445e-06
        },
        allocator,
    )
}
