use std::alloc::Allocator;

use super::{
    abs_max, double_to_sign, dummy_abs_max, predicates, ExplicitPoint3D, GenericNum,
    ImplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI, Orientation, Point3D,
};

/// Computes the orientation of the 3D points `pa`, `pb`, `pc` and `pd`.
/// If `pd` is above the plane defined by `pa`, `pb` and `pc`, then the
/// orientation is negative. If `pd` is below the plane, then the
/// orientation is positive. If `pd` is coplanar with `pa`, `pb` and `pc`,
/// then the orientation is zero.
pub fn orient3d<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    pd: &Point3D,
    bump: A,
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
            orient3d_llll(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lllt(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_llte(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lllt(pb, pa, pd, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pa, pb, pc, pd, bump)
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
            orient3d_lllt(pc, pd, pa, pb, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pc, pa, pb, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltte(pa, pb, pc, pd, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lltt(pa, pd, pb, pc, bump)
        }
        (Point3D::LPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pa, pb, pc, pd, bump)
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
            orient3d_lllt(pd, pc, pb, pa, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lltt(pb, pc, pa, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ltte(pb, pc, pa, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lltt(pb, pd, pc, pa, bump)
        }
        (Point3D::TPI(pa), Point3D::LPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pb, pa, pd, pc, bump)
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
            orient3d_lltt(pc, pd, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::LPI(pc), Point3D::TPI(pd)) => {
            orient3d_lttt(pc, pd, pa, pb, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::Explicit(pd)) => {
            orient3d_ttte(pc, pa, pb, pd, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::LPI(pd)) => {
            orient3d_lttt(pd, pc, pb, pa, bump)
        }
        (Point3D::TPI(pa), Point3D::TPI(pb), Point3D::TPI(pc), Point3D::TPI(pd)) => {
            orient3d_tttt(pa, pb, pc, pd, bump)
        }
    }
}

/// Implicit-Explicit-Explicit-Explicit orient3d
#[inline(always)]
fn orient3d_leee<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_teee<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_llee<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_ltee<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_ttee<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_llle<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_llte<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointTPI,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_ltte<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_ttte<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    pd: &ExplicitPoint3D,
    bump: A,
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
fn orient3d_llll<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    pd: &ImplicitPointLPI,
    bump: A,
) -> Orientation {
    orient3d_iiii(
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
            epsilon * 1.164303613521164e-07
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_lllt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointLPI,
    pd: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient3d_iiii(
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
        },
        bump,
    )
}
#[inline(always)]
fn orient3d_lltt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointLPI,
    pc: &ImplicitPointTPI,
    pd: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient3d_iiii(
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
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon * 0.001770733197190587
        },
        bump,
    )
}

#[inline(always)]
fn orient3d_lttt<A: Allocator + Copy>(
    pa: &ImplicitPointLPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    pd: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient3d_iiii(
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
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 0.01883943108077826
        },
        bump,
    )
}

#[inline(always)]
fn orient3d_tttt<A: Allocator + Copy>(
    pa: &ImplicitPointTPI,
    pb: &ImplicitPointTPI,
    pc: &ImplicitPointTPI,
    pd: &ImplicitPointTPI,
    bump: A,
) -> Orientation {
    orient3d_iiii(
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
            epsilon *= epsilon;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon *= max_var;
            epsilon * 0.1952243033447331
        },
        bump,
    )
}

fn orient3d_ieee_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    l1z: &T,
    d1: &T,
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
        abs_max(&[cx, cy, cz, ax_cx, ay_cy, az_cz, bx_cx, by_cy, bz_cz])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_ieee<IP: ImplicitPoint3D, F: FnOnce(f64) -> f64, A: Allocator + Copy>(
    pa: &IP,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(pa_static) = pa.static_filter() {
        let (det, max_var) = orient3d_ieee_impl::<true, _, _>(
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
        let (det, _) = orient3d_ieee_impl::<false, _, _>(
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
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(pa_exact) = pa.exact(bump) {
        let (det, _) = orient3d_ieee_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            [pb.data[0]].to_vec_in(bump).into(),
            [pb.data[1]].to_vec_in(bump).into(),
            [pb.data[2]].to_vec_in(bump).into(),
            [pc.data[0]].to_vec_in(bump).into(),
            [pc.data[1]].to_vec_in(bump).into(),
            [pc.data[2]].to_vec_in(bump).into(),
            [pd.data[0]].to_vec_in(bump).into(),
            [pd.data[1]].to_vec_in(bump).into(),
            [pd.data[2]].to_vec_in(bump).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiee_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    l1z: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    l2z: &T,
    d2: &T,
    p3x: T,
    p3y: T,
    p3z: T,
    p4x: T,
    p4y: T,
    p4z: T,
    abs_max: F,
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
        abs_max(&[p4x, p4y, p4z, p3p4x, p3p4y, p3p4z])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_iiee<
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    pa: &IP1,
    pb: &IP2,
    pc: &ExplicitPoint3D,
    pd: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some((pa_static, pb_static)) = pa.static_filter().zip(pb.static_filter()) {
        let (det, max_var) = orient3d_iiee_impl::<true, _, _>(
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
        let (det, _) = orient3d_iiee_impl::<false, _, _>(
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
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some((pa_exact, pb_exact)) = pa.exact(bump).zip(pb.exact(bump)) {
        let (det, _) = orient3d_iiee_impl::<false, _, _>(
            &pa_exact.x,
            &pa_exact.y,
            &pa_exact.z,
            &pa_exact.d,
            &pb_exact.x,
            &pb_exact.y,
            &pb_exact.z,
            &pb_exact.d,
            [pc.data[0]].to_vec_in(bump).into(),
            [pc.data[1]].to_vec_in(bump).into(),
            [pc.data[2]].to_vec_in(bump).into(),
            [pd.data[0]].to_vec_in(bump).into(),
            [pd.data[1]].to_vec_in(bump).into(),
            [pd.data[2]].to_vec_in(bump).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiie_impl<const NEED_MAX: bool, T: GenericNum, F: FnOnce(&[T]) -> Option<T>>(
    l1x: &T,
    l1y: &T,
    l1z: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    l2z: &T,
    d2: &T,
    l3x: &T,
    l3y: &T,
    l3z: &T,
    d3: &T,
    p4x: T,
    p4y: T,
    p4z: T,
    abs_max: F,
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
        abs_max(&[p4x, p4y, p4z])
    } else {
        None
    };
    (m012, max_var)
}

fn orient3d_iiie<
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    IP3: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    pa: &IP1,
    pb: &IP2,
    pc: &IP3,
    pd: &ExplicitPoint3D,
    static_filter_func: F,
    bump: A,
) -> Orientation {
    if let Some(((pa_static, pb_static), pc_static)) = pa
        .static_filter()
        .zip(pb.static_filter())
        .zip(pc.static_filter())
    {
        let (det, max_var) = orient3d_iiie_impl::<true, _, _>(
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
        let (det, _) = orient3d_iiie_impl::<false, _, _>(
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
        );
        if det.not_zero() {
            if det.positive() {
                return Orientation::Positive;
            } else if det.negative() {
                return Orientation::Negative;
            }
        }
    }

    if let Some(((pa_exact, pb_exact), pc_exact)) =
        pa.exact(bump).zip(pb.exact(bump)).zip(pc.exact(bump))
    {
        let (det, _) = orient3d_iiie_impl::<false, _, _>(
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
            [pd.data[0]].to_vec_in(bump).into(),
            [pd.data[1]].to_vec_in(bump).into(),
            [pd.data[2]].to_vec_in(bump).into(),
            dummy_abs_max,
        );
        return double_to_sign(*det.last().unwrap());
    }

    Orientation::Undefined
}

fn orient3d_iiii_impl<T: GenericNum>(
    l1x: &T,
    l1y: &T,
    l1z: &T,
    d1: &T,
    l2x: &T,
    l2y: &T,
    l2z: &T,
    d2: &T,
    l3x: &T,
    l3y: &T,
    l3z: &T,
    d3: &T,
    l4x: &T,
    l4y: &T,
    l4z: &T,
    d4: &T,
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
    IP1: ImplicitPoint3D,
    IP2: ImplicitPoint3D,
    IP3: ImplicitPoint3D,
    IP4: ImplicitPoint3D,
    F: FnOnce(f64) -> f64,
    A: Allocator + Copy,
>(
    pa: &IP1,
    pb: &IP2,
    pc: &IP3,
    pd: &IP4,
    static_filter_func: F,
    bump: A,
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

    if let Some((((pa_exact, pb_exact), pc_exact), pd_exact)) = pa
        .exact(bump)
        .zip(pb.exact(bump))
        .zip(pc.exact(bump))
        .zip(pd.exact(bump))
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
