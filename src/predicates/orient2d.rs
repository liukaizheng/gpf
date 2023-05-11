use bumpalo::{collections::Vec, vec, Bump};

use super::{
    abs_max, double_to_sign, dummy_abs_max, predicates, sign_reverse, ExplicitPoint3D, GenericNum,
    ImplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI, Orientation, Point3D,
};

/// Compute the orientation of the 3D points `pa`, `pb` and `pc`.
/// If `pc` is above the line defined by `pa` and `pb`, then the
/// orientation is negative. If `pc` is below the line, then the
/// orientation is positive. If `pc` is collinear with `pa` and `pb`,
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
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_xy(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_xy(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_xy(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_xy(pc, pa, pb),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_xy(pa, pb, pc),
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

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_xy(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_xy(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_xy(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_xy(pc, pa, pb),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_xy(pa, pb, pc),
    }
}

pub fn orient2d_zx<'a, 'b: 'a>(
    pa: &'a Point3D<'b>,
    pb: &'a Point3D<'b>,
    pc: &'a Point3D<'b>,
    bump: &'b Bump,
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
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_zx(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_zx(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_zx(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_zx(pc, pa, pb),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_zx(pa, pb, pc),
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

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_zx(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_zx(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_zx(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_zx(pc, pa, pb),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_zx(pa, pb, pc),
    }
}

pub fn orient2d_yz<'a, 'b: 'a>(
    pa: &'a Point3D<'b>,
    pb: &'a Point3D<'b>,
    pc: &'a Point3D<'b>,
    bump: &'b Bump,
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
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_lll_yz(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_llt_yz(pa, pb, pc),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_lte_yz(pa, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_llt_yz(pc, pa, pb),
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ltt_yz(pa, pb, pc),
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

        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc)) => orient2d_llt_yz(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc)) => orient2d_ltt_yz(pb, pc, pa),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc)) => {
            orient2d_tte_yz(pa, pb, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc)) => orient2d_ltt_yz(pc, pa, pb),
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc)) => orient2d_ttt_yz(pa, pb, pc),
    }
}

/// Implicit-Explicit-Explicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lee_xy<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lee::<'_, '_, 2>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_xy<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tee::<'_, '_, 2>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lle_xy<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointLPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lle::<'_, '_, 2>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_xy<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lte::<'_, '_, 2>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_xy<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tte::<'_, '_, 2>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for xy plane
#[inline(always)]
fn orient2d_lll_xy<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
) -> Orientation {
    orient2d_lll::<'_, '_, 2>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_llt_xy<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_llt::<'_, '_, 2>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ltt_xy<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ltt::<'_, '_, 2>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ttt_xy<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ttt::<'_, '_, 2>(pa, pb, pc)
}

/// Implicit-Explicit-Explicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lee_zx<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lee::<'_, '_, 1>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_zx<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tee::<'_, '_, 1>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lle_zx<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointLPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lle::<'_, '_, 1>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_zx<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lte::<'_, '_, 1>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_zx<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tte::<'_, '_, 1>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for zx plane
#[inline(always)]
fn orient2d_lll_zx<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
) -> Orientation {
    orient2d_lll::<'_, '_, 1>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_llt_zx<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_llt::<'_, '_, 1>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ltt_zx<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ltt::<'_, '_, 1>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ttt_zx<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ttt::<'_, '_, 1>(pa, pb, pc)
}

/// Implicit-Explicit-Explicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lee_yz<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lee::<'_, '_, 0>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tee_yz<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tee::<'_, '_, 0>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Explicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lle_yz<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointLPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lle::<'_, '_, 0>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_lte_yz<'b>(
    pa: &ImplicitPointLPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_lte::<'_, '_, 0>(pa, pb, pc, bump)
}
#[inline(always)]
fn orient2d_tte_yz<'b>(
    pa: &ImplicitPointTPI<'b>,
    pb: &ImplicitPointTPI<'b>,
    pc: &ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_tte::<'_, '_, 0>(pa, pb, pc, bump)
}

/// Implicit-Implicit-Implicit Orient2d for yz plane
#[inline(always)]
fn orient2d_lll_yz<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
) -> Orientation {
    orient2d_lll::<'_, '_, 0>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_llt_yz<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_llt::<'_, '_, 0>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ltt_yz<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ltt::<'_, '_, 0>(pa, pb, pc)
}
#[inline(always)]
fn orient2d_ttt_yz<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_ttt::<'_, '_, 0>(pa, pb, pc)
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
fn orient2d_iee<'a, 'b: 'a, const AXIS: u32, IP: ImplicitPoint3D<'b>, F: FnOnce(f64) -> f64>(
    p1: &'a IP,
    p2: &'a ExplicitPoint3D,
    p3: &'a ExplicitPoint3D,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation {
    if let Some(p1_static) = p1.static_filter() {
        let ret = if AXIS == 2 {
            orient2d_iee_impl::<'_, '_, true, _, _>(
                &p1_static.0.x,
                &p1_static.0.y,
                &p1_static.0.d,
                p2.data[0],
                p2.data[1],
                p3.data[0],
                p3.data[1],
                abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<'_, '_, true, _, _>(
                &p1_static.0.z,
                &p1_static.0.x,
                &p1_static.0.d,
                p2.data[2],
                p2.data[0],
                p3.data[2],
                p3.data[0],
                abs_max,
                bump,
            )
        } else {
            orient2d_iee_impl::<'_, '_, true, _, _>(
                &p1_static.0.y,
                &p1_static.0.z,
                &p1_static.0.d,
                p2.data[1],
                p2.data[2],
                p3.data[1],
                p3.data[2],
                abs_max,
                bump,
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
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_dynamic.x,
                &p1_dynamic.y,
                &p1_dynamic.d,
                p2.data[0].into(),
                p2.data[1].into(),
                p3.data[0].into(),
                p3.data[1].into(),
                dummy_abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_dynamic.z,
                &p1_dynamic.x,
                &p1_dynamic.d,
                p2.data[2].into(),
                p2.data[0].into(),
                p3.data[2].into(),
                p3.data[0].into(),
                dummy_abs_max,
                bump,
            )
        } else {
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_dynamic.y,
                &p1_dynamic.z,
                &p1_dynamic.d,
                p2.data[1].into(),
                p2.data[2].into(),
                p3.data[1].into(),
                p3.data[2].into(),
                dummy_abs_max,
                bump,
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

    if let Some(p1_exact) = p1.exact() {
        let (det, _) = if AXIS == 2 {
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_exact.x,
                &p1_exact.y,
                &p1_exact.d,
                vec![in bump; p2.data[0]].into(),
                vec![in bump; p2.data[1]].into(),
                vec![in bump; p3.data[0]].into(),
                vec![in bump; p3.data[1]].into(),
                dummy_abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_exact.z,
                &p1_exact.x,
                &p1_exact.d,
                vec![in bump; p2.data[2]].into(),
                vec![in bump; p2.data[0]].into(),
                vec![in bump; p3.data[2]].into(),
                vec![in bump; p3.data[0]].into(),
                dummy_abs_max,
                bump,
            )
        } else {
            orient2d_iee_impl::<'_, '_, false, _, _>(
                &p1_exact.y,
                &p1_exact.z,
                &p1_exact.d,
                vec![in bump; p2.data[1]].into(),
                vec![in bump; p2.data[2]].into(),
                vec![in bump; p3.data[1]].into(),
                vec![in bump; p3.data[2]].into(),
                dummy_abs_max,
                bump,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient2d_lee<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ExplicitPoint3D,
    p3: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iee::<'_, '_, AXIS, _, _>(
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

fn orient2d_tee<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointTPI<'b>,
    p2: &'a ExplicitPoint3D,
    p3: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iee::<'_, '_, AXIS, _, _>(
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

fn orient2d_iie_impl<
    'a,
    'b: 'a,
    const NEED_MAX: bool,
    T: 'b + GenericNum,
    F: FnOnce(Vec<'b, T>) -> Option<T>,
>(
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
        abs_max(vec![in bump; op3x, op3y])
    } else {
        None
    };
    (det, max_var)
}

fn orient2d_iie<
    'a,
    'b: 'a,
    const AXIS: u32,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    F: FnOnce(f64) -> f64,
>(
    p1: &'a IP1,
    p2: &'a IP2,
    p3: &'a ExplicitPoint3D,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation {
    if let Some((p1_static, p2_static)) = p1.static_filter().zip(p2.static_filter()) {
        let ret = if AXIS == 2 {
            orient2d_iie_impl::<'_, '_, true, _, _>(
                &p1_static.0.x,
                &p1_static.0.y,
                &p1_static.0.d,
                &p2_static.0.x,
                &p2_static.0.y,
                &p2_static.0.d,
                p3.data[0],
                p3.data[1],
                abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<'_, '_, true, _, _>(
                &p1_static.0.z,
                &p1_static.0.x,
                &p1_static.0.d,
                &p2_static.0.z,
                &p2_static.0.x,
                &p2_static.0.d,
                p3.data[2],
                p3.data[0],
                abs_max,
                bump,
            )
        } else {
            orient2d_iie_impl::<'_, '_, true, _, _>(
                &p1_static.0.y,
                &p1_static.0.z,
                &p1_static.0.d,
                &p2_static.0.y,
                &p2_static.0.z,
                &p2_static.0.d,
                p3.data[1],
                p3.data[2],
                abs_max,
                bump,
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
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_dynamic.x,
                &p1_dynamic.y,
                &p1_dynamic.d,
                &p2_dynamic.x,
                &p2_dynamic.y,
                &p2_dynamic.d,
                p3.data[0].into(),
                p3.data[1].into(),
                dummy_abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_dynamic.z,
                &p1_dynamic.x,
                &p1_dynamic.d,
                &p2_dynamic.z,
                &p2_dynamic.x,
                &p2_dynamic.d,
                p3.data[2].into(),
                p3.data[0].into(),
                dummy_abs_max,
                bump,
            )
        } else {
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_dynamic.y,
                &p1_dynamic.z,
                &p1_dynamic.d,
                &p2_dynamic.y,
                &p2_dynamic.z,
                &p2_dynamic.d,
                p3.data[1].into(),
                p3.data[2].into(),
                dummy_abs_max,
                bump,
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

    if let Some((p1_exact, p2_exact)) = p1.exact().zip(p2.exact()) {
        let (det, _) = if AXIS == 2 {
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_exact.x,
                &p1_exact.y,
                &p1_exact.d,
                &p2_exact.x,
                &p2_exact.y,
                &p2_exact.d,
                vec![in bump; p3.data[0]].into(),
                vec![in bump; p3.data[1]].into(),
                dummy_abs_max,
                bump,
            )
        } else if AXIS == 1 {
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_exact.z,
                &p1_exact.x,
                &p1_exact.d,
                &p2_exact.z,
                &p2_exact.x,
                &p2_exact.d,
                vec![in bump; p3.data[2]].into(),
                vec![in bump; p3.data[0]].into(),
                dummy_abs_max,
                bump,
            )
        } else {
            orient2d_iie_impl::<'_, '_, false, _, _>(
                &p1_exact.y,
                &p1_exact.z,
                &p1_exact.d,
                &p2_exact.y,
                &p2_exact.z,
                &p2_exact.d,
                vec![in bump; p3.data[1]].into(),
                vec![in bump; p3.data[2]].into(),
                dummy_abs_max,
                bump,
            )
        };
        return double_to_sign(*det.last().unwrap());
    }
    Orientation::Undefined
}

#[inline(always)]
fn orient2d_lle<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointLPI<'b>,
    p3: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie::<'_, '_, AXIS, _, _, _>(
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
fn orient2d_lte<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie::<'_, '_, AXIS, _, _, _>(
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
fn orient2d_tte<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointTPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient2d_iie::<'_, '_, AXIS, _, _, _>(
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
    'a,
    'b: 'a,
    const AXIS: u32,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    IP3: ImplicitPoint3D<'b>,
    F: FnOnce(f64) -> f64,
>(
    p1: &'a IP1,
    p2: &'a IP2,
    p3: &'a IP3,
    static_filter_func: F,
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

    if let Some(((p1_exact, p2_exact), p3_exact)) = p1.exact().zip(p2.exact()).zip(p3.exact()) {
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
fn orient2d_lll<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointLPI<'b>,
    p3: &'a ImplicitPointLPI<'b>,
) -> Orientation {
    orient2d_iii::<'_, '_, AXIS, _, _, _, _>(p1, p2, p3, |max_var| {
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
fn orient2d_llt<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointLPI<'b>,
    p3: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii::<'_, '_, AXIS, _, _, _, _>(p1, p2, p3, |max_var| {
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
fn orient2d_ltt<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointLPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii::<'_, '_, AXIS, _, _, _, _>(p1, p2, p3, |max_var| {
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
fn orient2d_ttt<'a, 'b: 'a, const AXIS: u32>(
    p1: &'a ImplicitPointTPI<'b>,
    p2: &'a ImplicitPointTPI<'b>,
    p3: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient2d_iii::<'_, '_, AXIS, _, _, _, _>(p1, p2, p3, |max_var| {
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
