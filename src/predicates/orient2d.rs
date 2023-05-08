use bumpalo::Bump;

use super::{double_to_sign, predicates, ExplicitPoint3D, Orientation, Point3D, ImplicitPointLPI};

/// Computes the orientation of the 3D points `pa`, `pb`, `pc` and `pd`.
/// If `pd` is above the plane defined by `pa`, `pb` and `pc`, then the
/// orientation is negative. If `pd` is below the plane, then the
/// orientation is positive. If `pd` is coplanar with `pa`, `pb` and `pc`,
/// then the orientation is zero.
pub fn orient2d(pa: &Point3D, pb: &Point3D, pc: &Point3D, bump: &Bump) -> Orientation {
    match (pa, pb, pc) {
        (Point3D::Explicit(pa), Point3D::Explicit(pb), Point3D::Explicit(pc)) => {
            double_to_sign(predicates::orient2d(pa, pb, pc, bump))
        }
        (Point3D::Explicit(_), Point3D::Explicit(_), Point3D::LPI(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::Explicit(_), Point3D::TPI(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::LPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::LPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::LPI(_), Point3D::TPI(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::TPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::TPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::Explicit(_), Point3D::TPI(_), Point3D::TPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::Explicit(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::LPI(_), Point3D::Explicit(_), Point3D::LPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::Explicit(_), Point3D::TPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::LPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::LPI(_), Point3D::LPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::LPI(_), Point3D::TPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::TPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::LPI(_), Point3D::TPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::LPI(_), Point3D::TPI(_), Point3D::TPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::Explicit(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::TPI(_), Point3D::Explicit(_), Point3D::LPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::Explicit(_), Point3D::TPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::LPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::TPI(_), Point3D::LPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::LPI(_), Point3D::TPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::TPI(_), Point3D::Explicit(_)) => todo!(),
        (Point3D::TPI(_), Point3D::TPI(_), Point3D::LPI(_)) => todo!(),
        (Point3D::TPI(_), Point3D::TPI(_), Point3D::TPI(_)) => todo!(),
    }
}

fn orient2d_lee(pa: ImplicitPointLPI, pb: &[f64], pc: &[f64]) -> Orientation {
	Orientation::Zero
}
