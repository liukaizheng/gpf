use bumpalo::{collections::Vec, vec, Bump};

use super::{
    abs_max, double_to_sign, dummy_abs_max, predicates, ExplicitPoint3D, GenericNum,
    ImplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI, Orientation, Point3D,
};

/// Computes the orientation of the 3D points `pa`, `pb`, `pc` and `pd`.
/// If `pd` is above the plane defined by `pa`, `pb` and `pc`, then the
/// orientation is negative. If `pd` is below the plane, then the
/// orientation is positive. If `pd` is coplanar with `pa`, `pb` and `pc`,
/// then the orientation is zero.
pub fn orient3d<'a, 'b: 'a>(
    pa: &'a Point3D<'b>,
    pb: &'a Point3D<'b>,
    pc: &'a Point3D<'b>,
    pd: &'a Point3D<'b>,
    bump: &'b Bump,
) -> Orientation {
    match (pa, pb, pc, pd) {
        (
            Point3D::Explicit(pa),
            Point3D::Explicit(pb),
            Point3D::Explicit(pc),
            Point3D::Explicit(pd),
        ) => double_to_sign(-predicates::orient3d(pa, pb, pc, pd, bump)),
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_leee(pd, pa, pc, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_teee(pd, pa, pc, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_leee(pc, pd, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llee(pc, pd, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_ltee(pc, pd, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_teee(pc, pd, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_ltee(pd, pc, pb, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_ttee(pc, pd, pa, pb, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_leee(pb, pc, pa, pd, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_llee(pb, pd, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ltee(pb, pd, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llee(pb, pc, pa, pd, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llle(pb, pd, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_llte(pc, pb, pd, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pb, pc, pa, pd, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pb, pd, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pb, pd, pc, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_teee(pb, pc, pa, pd, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_ltee(pd, pb, pa, pc, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ttee(pd, pb, pa, pc, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pc, pb, pd, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pd, pc, pb, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pc, pb, pd, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ttee(pb, pc, pa, pd, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_ltte(pd, pc, pb, pa, bump)
        }
        (Point3D::Explicit(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_ttte(pd, pc, pb, pa, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_leee(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_llee(pd, pa, pc, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ltee(pa, pd, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llee(pa, pc, pd, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llle(pa, pc, pd, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_llte(pa, pc, pd, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pa, pc, pd, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pd, pa, pc, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pa, pc, pd, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_llee(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_llle(pa, pd, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_llte(pb, pa, pd, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llle(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llll(pa, pb, pc, pd)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lllt(pa, pb, pc, pd)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llte(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lllt(pb, pa, pd, pc)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pa, pb, pc, pd)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pa, pd, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pa, pd, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llte(pc, pa, pb, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_lllt(pc, pd, pa, pb)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pc, pa, pb, pd)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltte(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lltt(pa, pd, pb, pc)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pa, pb, pc, pd)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_teee(pa, pb, pc, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_ltee(pd, pa, pc, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ttee(pd, pa, pc, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pc, pa, pb, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pc, pd, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pc, pd, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ttee(pc, pa, pb, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_ltte(pd, pa, pc, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::Explicit(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_ttte(pd, pa, pc, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_ltee(pb, pa, pd, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_llte(pd, pb, pa, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ltte(pb, pa, pd, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llte(pb, pc, pa, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_lllt(pd, pc, pb, pa)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pb, pc, pa, pd)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltte(pb, pc, pa, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lltt(pb, pd, pc, pa)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pb, pa, pd, pc)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::Explicit(pd)) => {
            orient3d_ttee(pb, pa, pd, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::LPI(pd)) => {
            orient3d_ltte(pd, pb, pa, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::Explicit(pc), Point3D::TPI(pd)) => {
            orient3d_ttte(pd, pb, pa, pc, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltte(pc, pa, pb, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::LPI(pd)) => {
            orient3d_lltt(pc, pd, pa, pb)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pc, pd, pa, pb)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ttte(pc, pa, pb, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lttt(pd, pc, pb, pa)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_tttt(pa, pb, pc, pd)
        }
    }
}

/// Implicit-Explicit-Explicit-Explicit orient3d
#[inline(always)]
fn orient3d_leee<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ExplicitPoint3D,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_ieee(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 1.861039534284405e-13
        },
        bump,
    )
}

#[inline(always)]
fn orient3d_teee<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ExplicitPoint3D,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_ieee(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon * 3.070283610684406e-12
        },
        bump,
    )
}

/// Implicit-Implicit-Explicit-Explicit orient3d
#[inline(always)]
fn orient3d_llee<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiee(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon * 5.12855469897434e-12
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_ltee<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiee(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 7.437036403379365e-11
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_ttee<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiee(
        pa,
        pb,
        pc,
        pd,
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
            epsilon *= max_var;
            epsilon * 1.036198238324465e-09
        },
        bump,
    )
}

/// Implicit-Implicit-Implicit-Explicit orient3d
#[inline(always)]
fn orient3d_llle<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiie(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 1.270161397934349e-10
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_llte<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiie(
        pa,
        pb,
        pc,
        pd,
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
            epsilon *= max_var;
            epsilon * 1.7060943907632e-09
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_ltte<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiie(
        pa,
        pb,
        pc,
        pd,
        |max_var| {
            let mut epsilon = max_var;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 2.211968919141342e-08
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_ttte<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ExplicitPoint3D,
    bump: &'b Bump,
) -> Orientation {
    orient3d_iiie(
        pa,
        pb,
        pc,
        pd,
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
            epsilon * 2.808754828720361e-07
        },
        bump,
    )
}

/// Implicit-Implicit-Implicit-Explicit orient3d
#[inline(always)]
fn orient3d_llll<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
    pd: &'a ImplicitPointLPI<'b>,
) -> Orientation {
    orient3d_iiii(pa, pb, pc, pd, |max_var| {
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
        epsilon * 1.164303613521164e-07
    })
}
#[inline(always)]
fn orient3d_lllt<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointLPI<'b>,
    pd: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient3d_iiii(pa, pb, pc, pd, |max_var| {
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
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon * 0.0001675978376241023
    })
}
#[inline(always)]
fn orient3d_lltt<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointLPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient3d_iiii(pa, pb, pc, pd, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon * 0.001770733197190587
    })
}

#[inline(always)]
fn orient3d_lttt<'a, 'b: 'a>(
    pa: &'a ImplicitPointLPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient3d_iiii(pa, pb, pc, pd, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon *= max_var;
        epsilon * 0.01883943108077826
    })
}

#[inline(always)]
fn orient3d_tttt<'a, 'b: 'a>(
    pa: &'a ImplicitPointTPI<'b>,
    pb: &'a ImplicitPointTPI<'b>,
    pc: &'a ImplicitPointTPI<'b>,
    pd: &'a ImplicitPointTPI<'b>,
) -> Orientation {
    orient3d_iiii(pa, pb, pc, pd, |max_var| {
        let mut epsilon = max_var;
        epsilon *= epsilon;
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
        epsilon * 0.1952243033447331
    })
}

fn orient3d_ieee_impl<
    'a,
    'b: 'a,
    const NEED_MAX: bool,
    T: 'b + GenericNum,
    F: FnOnce(Vec<'b, T>) -> Option<T>,
>(
    l1x: &'a T,
    l1y: &'a T,
    l1z: &'a T,
    d1: &'a T,
    ax: T,
    ay: T,
    az: T,
    bx: T,
    by: T,
    bz: T,
    cx: T,
    cy: T,
    cz: T,
    abs_max: F,
    bump: &'b Bump,
) -> (T, Option<T>) {
    let dcx = d1 * &cx;
    let dcy = d1 * &cy;
    let dcz = d1 * &cz;
    let ix_cx = l1x - dcx;
    let iy_cy = l1y - dcy;
    let ax_cx = ax - &cx;
    let ay_cy = ay - &cy;
    let az_cz = az - &cz;
    let iz_cz = l1z - dcz;
    let bx_cx = bx - &cx;
    let by_cy = by - &cy;
    let bz_cz = bz - &cz;
    let tmc_a = &ix_cx * &ay_cy;
    let tmc_b = &iy_cy * &ax_cx;
    let m01 = tmc_a - tmc_b;
    let tmi_a = ix_cx * &az_cz;
    let tmi_b = &iz_cz * &ax_cx;
    let m02 = tmi_a - tmi_b;
    let tma_a = iy_cy * &az_cz;
    let tma_b = iz_cz * &ay_cy;
    let m12 = tma_a - tma_b;
    let mt1 = m01 * &bz_cz;
    let mt2 = m02 * &by_cy;
    let mt3 = m12 * &bx_cx;
    let mtt = mt2 - mt1;
    let m012 = mtt - mt3;

    let max_var = if NEED_MAX {
        abs_max(vec![
            in bump; cx, cy, cz,
            ax_cx, ay_cy, az_cz,
            bx_cx, by_cy, bz_cz,
        ])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_ieee<'a, 'b: 'a, IP: ImplicitPoint3D<'b>, F: FnOnce(f64) -> f64>(
    pa: &'a IP,
    pb: &'a ExplicitPoint3D,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation {
    if let Some(pa_static) = pa.static_filter() {
        let (det, max_var) = orient3d_ieee_impl::<'_, '_, true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.z,
            &pa_static.0.d,
            pb.data[0],
            pb.data[1],
            pb.data[2],
            pc.data[0],
            pc.data[1],
            pc.data[2],
            pd.data[0],
            pd.data[1],
            pd.data[2],
            abs_max,
            bump,
        );
        let max_var = max_var.unwrap().max(pa_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(pa_dynamic) = pa.dynamic_filter() {
        let (det, _) = orient3d_ieee_impl::<'_, '_, false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.z,
            &pa_dynamic.d,
            pb.data[0].into(),
            pb.data[1].into(),
            pb.data[2].into(),
            pc.data[0].into(),
            pc.data[1].into(),
            pc.data[2].into(),
            pd.data[0].into(),
            pd.data[1].into(),
            pd.data[2].into(),
            dummy_abs_max,
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

    if let Some(pa_exact) = pa.exact() {
        let (det, _) = orient3d_ieee_impl::<'_, '_, false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            vec![in bump; pb.data[0]].into(),
            vec![in bump; pb.data[1]].into(),
            vec![in bump; pb.data[2]].into(),
            vec![in bump; pc.data[0]].into(),
            vec![in bump; pc.data[1]].into(),
            vec![in bump; pc.data[2]].into(),
            vec![in bump; pd.data[0]].into(),
            vec![in bump; pd.data[1]].into(),
            vec![in bump; pd.data[2]].into(),
            dummy_abs_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiee_impl<
    'a,
    'b: 'a,
    const NEED_MAX: bool,
    T: 'b + GenericNum,
    F: FnOnce(Vec<'b, T>) -> Option<T>,
>(
    l1x: &'a T,
    l1y: &'a T,
    l1z: &'a T,
    d1: &'a T,
    l2x: &'a T,
    l2y: &'a T,
    l2z: &'a T,
    d2: &'a T,
    p3x: T,
    p3y: T,
    p3z: T,
    p4x: T,
    p4y: T,
    p4z: T,
    abs_max: F,
    bump: &'b Bump,
) -> (T, Option<T>) {
    let d1p4x = d1 * &p4x;
    let d1p4y = d1 * &p4y;
    let d1p4z = d1 * &p4z;
    let d2p4x = d2 * &p4x;
    let d2p4y = d2 * &p4y;
    let d2p4z = d2 * &p4z;
    let p1p4x = l1x - d1p4x;
    let p1p4y = l1y - d1p4y;
    let p1p4z = l1z - d1p4z;
    let p2p4x = l2x - d2p4x;
    let p2p4y = l2y - d2p4y;
    let p2p4z = l2z - d2p4z;
    let p3p4x = p3x - &p4x;
    let p3p4y = p3y - &p4y;
    let p3p4z = p3z - &p4z;
    let tmc_a = &p1p4x * &p2p4y;
    let tmc_b = &p1p4y * &p2p4x;
    let m01 = tmc_a - tmc_b;
    let tmi_a = p1p4x * &p2p4z;
    let tmi_b = &p1p4z * &p2p4x;
    let m02 = tmi_a - tmi_b;
    let tma_a = p1p4y * p2p4z;
    let tma_b = p1p4z * p2p4y;
    let m12 = tma_a - tma_b;
    let mt1 = m01 * &p3p4z;
    let mt2 = m02 * &p3p4y;
    let mt3 = m12 * &p3p4x;
    let mtt = mt2 - mt1;
    let m012 = mtt - mt3;

    let max_var = if NEED_MAX {
        abs_max(vec![in bump; p4x, p4y, p4z, p3p4x, p3p4y, p3p4z])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_iiee<
    'a,
    'b: 'a,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    F: FnOnce(f64) -> f64,
>(
    pa: &'a IP1,
    pb: &'a IP2,
    pc: &'a ExplicitPoint3D,
    pd: &'a ExplicitPoint3D,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation {
    if let Some((pa_static, pb_static)) = pa.static_filter().zip(pb.static_filter()) {
        let (det, max_var) = orient3d_iiee_impl::<'_, '_, true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.z,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.z,
            &pb_static.0.d,
            pc.data[0],
            pc.data[1],
            pc.data[2],
            pd.data[0],
            pd.data[1],
            pd.data[2],
            abs_max,
            bump,
        );
        let max_var = max_var.unwrap().max(pa_static.1).max(pb_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some((pa_dynamic, pb_dynamic)) = pa.dynamic_filter().zip(pb.dynamic_filter()) {
        let (det, _) = orient3d_iiee_impl::<'_, '_, false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.z,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.z,
            &pb_dynamic.d,
            pc.data[0].into(),
            pc.data[1].into(),
            pc.data[2].into(),
            pd.data[0].into(),
            pd.data[1].into(),
            pd.data[2].into(),
            dummy_abs_max,
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

    if let Some((pa_exact, pb_exact)) = pa.exact().zip(pb.exact()) {
        let (det, _) = orient3d_iiee_impl::<'_, '_, false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.z,
            &pb_exact.d,
            vec![in bump; pc.data[0]].into(),
            vec![in bump; pc.data[1]].into(),
            vec![in bump; pc.data[2]].into(),
            vec![in bump; pd.data[0]].into(),
            vec![in bump; pd.data[1]].into(),
            vec![in bump; pd.data[2]].into(),
            dummy_abs_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiie_impl<
    'a,
    'b: 'a,
    const NEED_MAX: bool,
    T: 'b + GenericNum,
    F: FnOnce(Vec<'b, T>) -> Option<T>,
>(
    l1x: &'a T,
    l1y: &'a T,
    l1z: &'a T,
    d1: &'a T,
    l2x: &'a T,
    l2y: &'a T,
    l2z: &'a T,
    d2: &'a T,
    l3x: &'a T,
    l3y: &'a T,
    l3z: &'a T,
    d3: &'a T,
    p4x: T,
    p4y: T,
    p4z: T,
    abs_max: F,
    bump: &'b Bump,
) -> (T, Option<T>) {
    let d1p4x = d1 * &p4x;
    let d1p4y = d1 * &p4y;
    let d1p4z = d1 * &p4z;
    let d2p4x = d2 * &p4x;
    let d2p4y = d2 * &p4y;
    let d2p4z = d2 * &p4z;
    let d3p4x = d3 * &p4x;
    let d3p4y = d3 * &p4y;
    let d3p4z = d3 * &p4z;
    let p1p4x = l1x - d1p4x;
    let p1p4y = l1y - d1p4y;
    let p1p4z = l1z - d1p4z;
    let p2p4x = l2x - d2p4x;
    let p2p4y = l2y - d2p4y;
    let p2p4z = l2z - d2p4z;
    let p3p4x = l3x - d3p4x;
    let p3p4y = l3y - d3p4y;
    let p3p4z = l3z - d3p4z;
    let tmc_a = &p1p4x * &p2p4y;
    let tmc_b = &p1p4y * &p2p4x;
    let m01 = tmc_a - tmc_b;
    let tmi_a = p1p4x * &p2p4z;
    let tmi_b = &p1p4z * p2p4x;
    let m02 = tmi_a - tmi_b;
    let tma_a = p1p4y * p2p4z;
    let tma_b = p1p4z * p2p4y;
    let m12 = tma_a - tma_b;
    let mt1 = m01 * p3p4z;
    let mt2 = m02 * p3p4y;
    let mt3 = m12 * p3p4x;
    let mtt = mt2 - mt1;
    let m012 = mtt - mt3;
    let max_var = if NEED_MAX {
        abs_max(vec![in bump; p4x, p4y, p4z])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_iiie<
    'a,
    'b: 'a,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    IP3: ImplicitPoint3D<'b>,
    F: FnOnce(f64) -> f64,
>(
    pa: &'a IP1,
    pb: &'a IP2,
    pc: &'a IP3,
    pd: &'a ExplicitPoint3D,
    static_filter_func: F,
    bump: &'b Bump,
) -> Orientation {
    if let Some(((pa_static, pb_static), pc_static)) = pa
        .static_filter()
        .zip(pb.static_filter())
        .zip(pc.static_filter())
    {
        let (det, max_var) = orient3d_iiie_impl::<'_, '_, true, _, _>(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.z,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.z,
            &pb_static.0.d,
            &pc_static.0.x,
            &pc_static.0.y,
            &pc_static.0.z,
            &pc_static.0.d,
            pd.data[0],
            pd.data[1],
            pd.data[2],
            abs_max,
            bump,
        );
        let max_var = max_var
            .unwrap()
            .max(pa_static.1)
            .max(pb_static.1)
            .max(pc_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some(((pa_dynamic, pb_dynamic), pc_dynamic)) = pa
        .dynamic_filter()
        .zip(pb.dynamic_filter())
        .zip(pc.dynamic_filter())
    {
        let (det, _) = orient3d_iiie_impl::<'_, '_, false, _, _>(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.z,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.z,
            &pb_dynamic.d,
            &pc_dynamic.x,
            &pc_dynamic.y,
            &pc_dynamic.z,
            &pc_dynamic.d,
            pd.data[0].into(),
            pd.data[1].into(),
            pd.data[2].into(),
            dummy_abs_max,
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

    if let Some(((pa_exact, pb_exact), pc_exact)) = pa.exact().zip(pb.exact()).zip(pc.exact()) {
        let (det, _) = orient3d_iiie_impl::<'_, '_, false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.z,
            &pb_exact.d,
            &pc_exact.x,
            &pc_exact.y,
            &pc_exact.z,
            &pc_exact.d,
            vec![in bump; pd.data[0]].into(),
            vec![in bump; pd.data[1]].into(),
            vec![in bump; pd.data[2]].into(),
            dummy_abs_max,
            bump,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiii_impl<'a, T: GenericNum>(
    l1x: &'a T,
    l1y: &'a T,
    l1z: &'a T,
    d1: &'a T,
    l2x: &'a T,
    l2y: &'a T,
    l2z: &'a T,
    d2: &'a T,
    l3x: &'a T,
    l3y: &'a T,
    l3z: &'a T,
    d3: &'a T,
    l4x: &'a T,
    l4y: &'a T,
    l4z: &'a T,
    d4: &'a T,
) -> T {
    let d1p4x = d1 * l4x;
    let d1p4y = d1 * l4y;
    let d1p4z = d1 * l4z;
    let d2p4x = d2 * l4x;
    let d2p4y = d2 * l4y;
    let d2p4z = d2 * l4z;
    let d3p4x = d3 * l4x;
    let d3p4y = d3 * l4y;
    let d3p4z = d3 * l4z;
    let d4l1x = d4 * l1x;
    let d4l1y = d4 * l1y;
    let d4l1z = d4 * l1z;
    let d4l2x = d4 * l2x;
    let d4l2y = d4 * l2y;
    let d4l2z = d4 * l2z;
    let d4l3x = d4 * l3x;
    let d4l3y = d4 * l3y;
    let d4l3z = d4 * l3z;
    let p1p4x = d4l1x - d1p4x;
    let p1p4y = d4l1y - d1p4y;
    let p1p4z = d4l1z - d1p4z;
    let p2p4x = d4l2x - d2p4x;
    let p2p4y = d4l2y - d2p4y;
    let p2p4z = d4l2z - d2p4z;
    let p3p4x = d4l3x - d3p4x;
    let p3p4y = d4l3y - d3p4y;
    let p3p4z = d4l3z - d3p4z;
    let tmc_a = &p1p4x * &p2p4y;
    let tmc_b = &p1p4y * &p2p4x;
    let m01 = tmc_a - tmc_b;
    let tmi_a = p1p4x * &p2p4z;
    let tmi_b = &p1p4z * p2p4x;
    let m02 = tmi_a - tmi_b;
    let tma_a = p1p4y * p2p4z;
    let tma_b = p1p4z * p2p4y;
    let m12 = tma_a - tma_b;
    let mt1 = m01 * p3p4z;
    let mt2 = m02 * p3p4y;
    let mt3 = m12 * p3p4x;
    let mtt = mt2 - mt1;
    mtt - mt3
}

fn orient3d_iiii<
    'a,
    'b: 'a,
    IP1: ImplicitPoint3D<'b>,
    IP2: ImplicitPoint3D<'b>,
    IP3: ImplicitPoint3D<'b>,
    IP4: ImplicitPoint3D<'b>,
    F: FnOnce(f64) -> f64,
>(
    pa: &'a IP1,
    pb: &'a IP2,
    pc: &'a IP3,
    pd: &'a IP4,
    static_filter_func: F,
) -> Orientation {
    if let Some((((pa_static, pb_static), pc_static), pd_static)) = pa
        .static_filter()
        .zip(pb.static_filter())
        .zip(pc.static_filter())
        .zip(pd.static_filter())
    {
        let det = orient3d_iiii_impl(
            &pa_static.0.x,
            &pa_static.0.y,
            &pa_static.0.z,
            &pa_static.0.d,
            &pb_static.0.x,
            &pb_static.0.y,
            &pb_static.0.z,
            &pb_static.0.d,
            &pc_static.0.x,
            &pc_static.0.y,
            &pc_static.0.z,
            &pc_static.0.d,
            &pd_static.0.x,
            &pd_static.0.y,
            &pd_static.0.z,
            &pd_static.0.d,
        );
        let max_var = pa_static
            .1
            .max(pb_static.1)
            .max(pc_static.1)
            .max(pd_static.1);
        let epsilon = static_filter_func(max_var);
        if det > epsilon {
            return Orientation::Positive;
        } else if det < -epsilon {
            return Orientation::Negative;
        }
    }

    if let Some((((pa_dynamic, pb_dynamic), pc_dynamic), pd_dynamic)) = pa
        .dynamic_filter()
        .zip(pb.dynamic_filter())
        .zip(pc.dynamic_filter())
        .zip(pd.dynamic_filter())
    {
        let det = orient3d_iiii_impl(
            &pa_dynamic.x,
            &pa_dynamic.y,
            &pa_dynamic.z,
            &pa_dynamic.d,
            &pb_dynamic.x,
            &pb_dynamic.y,
            &pb_dynamic.z,
            &pb_dynamic.d,
            &pc_dynamic.x,
            &pc_dynamic.y,
            &pc_dynamic.z,
            &pc_dynamic.d,
            &pd_dynamic.x,
            &pd_dynamic.y,
            &pd_dynamic.z,
            &pd_dynamic.d,
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some((((pa_exact, pb_exact), pc_exact), pd_exact)) =
        pa.exact().zip(pb.exact()).zip(pc.exact()).zip(pd.exact())
    {
        let det = orient3d_iiii_impl(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.z,
            &pb_exact.d,
            &pc_exact.x,
            &pc_exact.y,
            &pc_exact.z,
            &pc_exact.d,
            &pd_exact.x,
            &pd_exact.y,
            &pd_exact.z,
            &pd_exact.d,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}
