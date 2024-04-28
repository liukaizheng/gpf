mod expansion_number;
mod generic_point;
mod interval_number;
pub mod orient2d;
pub mod orient3d;
mod predicates;
mod less_than;

use std::{
    alloc::Allocator, ops::{Add, Mul, Sub}
};

pub use expansion_number::*;
pub use generic_point::*;
pub use interval_number::*;
pub use orient2d::*;
pub use predicates::*;
pub use less_than::*;


#[derive(PartialEq, Eq, Clone, Copy, Hash, Debug)]
pub enum Orientation {
    Positive,
    Negative,
    Zero,
    Undefined,
}

trait GenericNum = Sized + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self>
where
    for<'a> Self:
        Add<&'a Self, Output = Self> + Sub<&'a Self, Output = Self> + Mul<&'a Self, Output = Self>,
    for<'a> &'a Self: Add<Output = Self>
        + Add<Self, Output = Self>
        + Sub<Output = Self>
        + Sub<Self, Output = Self>
        + Mul<Output = Self>
        + Mul<Self, Output = Self>;

#[inline(always)]
fn abs_max(vec: &[f64]) -> Option<f64> {
    vec.into_iter().map(|e| e.abs()).reduce(|acc, e| acc.max(e))
}

#[inline(always)]
fn dummy_abs_max<T>(_: &[T]) -> Option<T> {
    None
}

#[inline(always)]
pub fn double_to_sign(x: f64) -> Orientation {
    if x > 0.0 {
        Orientation::Positive
    } else if x < 0.0 {
        Orientation::Negative
    } else {
        Orientation::Zero
    }
}

#[inline(always)]
pub fn sign_reverse(ori: Orientation) -> Orientation {
    match ori {
        Orientation::Positive => Orientation::Negative,
        Orientation::Negative => Orientation::Positive,
        ori => ori,
    }
}

#[inline(always)]
pub fn sign_reversed(ori1: Orientation, ori2: Orientation) -> bool {
    (ori1 == Orientation::Positive && ori2 == Orientation::Negative)
        || (ori2 == Orientation::Positive && ori1 == Orientation::Negative)
}

#[inline(always)]
pub fn get_exponent(x: f64) -> i32 {
    if x == 0.0 {
        0
    } else {
        x.abs().to_bits().wrapping_shr(52) as i32 - 1023
    }
}

fn max_comp_in_tri_normal_impl<
    const NEED_MAX: bool,
    T: GenericNum,
    F: FnOnce(&[T]) -> Option<T>,
>(
    ov1x: T,
    ov1y: T,
    ov1z: T,
    ov2x: T,
    ov2y: T,
    ov2z: T,
    ov3x: T,
    ov3y: T,
    ov3z: T,
    abs_max: F,
) -> (T, T, T, Option<T>) {
    let v3x = ov3x - &ov2x;
    let v3y = ov3y - &ov2y;
    let v3z = ov3z - &ov2z;
    let v2x = ov2x - ov1x;
    let v2y = ov2y - ov1y;
    let v2z = ov2z - ov1z;
    let nvx1 = &v2y * &v3z;
    let nvx2 = &v2z * &v3y;
    let nvx = nvx1 - nvx2;
    let nvy1 = &v3x * &v2z;
    let nvy2 = &v3z * &v2x;
    let nvy = nvy1 - nvy2;
    let nvz1 = &v2x * &v3y;
    let nvz2 = &v2y * &v3x;
    let nvz = nvz1 - nvz2;
    let max_var = if NEED_MAX {
        abs_max(&[v3x, v3y, v3z, v2x, v2y, v2z])
    } else {
        None
    };
    (nvx, nvy, nvz, max_var)
}

#[inline(always)]
fn max_comp_in_tri_normal_filter(
    ov1x: f64,
    ov1y: f64,
    ov1z: f64,
    ov2x: f64,
    ov2y: f64,
    ov2z: f64,
    ov3x: f64,
    ov3y: f64,
    ov3z: f64,
) -> usize {
    let (nvx, nvy, nvz, max_var) = max_comp_in_tri_normal_impl::<true, _, _>(
        ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z, abs_max,
    );
    let max_var = max_var.unwrap();
    let epsilon = 8.88395e-016 * max_var * max_var;
    let (max_comp, max_val) = [nvx.abs(), nvy.abs(), nvz.abs()]
        .into_iter()
        .enumerate()
        .max_by(|(_, v1), (_, v2)| v1.total_cmp(v2))
        .unwrap();
    if max_val > epsilon {
        max_comp
    } else {
        3
    }
}

#[inline(always)]
fn max_comp_in_tri_normal_exact<A: Allocator + Copy>(
    ov1x: ExpansionNum<A>,
    ov1y: ExpansionNum<A>,
    ov1z: ExpansionNum<A>,
    ov2x: ExpansionNum<A>,
    ov2y: ExpansionNum<A>,
    ov2z: ExpansionNum<A>,
    ov3x: ExpansionNum<A>,
    ov3y: ExpansionNum<A>,
    ov3z: ExpansionNum<A>,
) -> usize {
    let (nvx, nvy, nvz, _) = max_comp_in_tri_normal_impl::<false, _, _>(
        ov1x,
        ov1y,
        ov1z,
        ov2x,
        ov2y,
        ov2z,
        ov3x,
        ov3y,
        ov3z,
        dummy_abs_max,
    );
    let (max_comp, _) = [nvx, nvy, nvz]
        .into_iter()
        .map(|x| x.vec.last().unwrap().abs())
        .enumerate()
        .max_by(|(_, v1), (_, v2)| v1.total_cmp(v2))
        .unwrap();
    max_comp
}

#[inline(always)]
pub fn max_comp_in_tri_normal<A: Allocator + Copy>(
    ov1: &[f64],
    ov2: &[f64],
    ov3: &[f64],
    allocator: A,
) -> usize {
    let max_comp = max_comp_in_tri_normal_filter(
        ov1[0], ov1[1], ov1[2], ov2[0], ov2[1], ov2[2], ov3[0], ov3[1], ov3[2],
    );
    if max_comp < 3 {
        max_comp
    } else {
        max_comp_in_tri_normal_exact(
            [ov1[0]].to_vec_in(allocator).into(),
            [ov1[1]].to_vec_in(allocator).into(),
            [ov1[2]].to_vec_in(allocator).into(),
            [ov2[0]].to_vec_in(allocator).into(),
            [ov2[1]].to_vec_in(allocator).into(),
            [ov2[2]].to_vec_in(allocator).into(),
            [ov3[0]].to_vec_in(allocator).into(),
            [ov3[1]].to_vec_in(allocator).into(),
            [ov3[2]].to_vec_in(allocator).into(),
        )
    }
}

pub fn point_in_inner_triangle<A: Allocator + Copy>(
    p: &[f64],
    v1: &[f64],
    v2: &[f64],
    v3: &[f64],
    allocator: A,
) -> bool {
    // Projection on (x,y)-plane -> p VS v1
    let mut o1 = double_to_sign(predicates::orient2d(p, v2, v3, allocator));
    let mut o2 = double_to_sign(predicates::orient2d(v1, v2, v3, allocator));
    let oo2 = o2;
    if o1 != o2 {
        return false;
    };

    // Projection on (y,z)-plane -> p VS v1
    o1 = double_to_sign(predicates::orient2d(&p[1..], &v2[1..], &v3[1..], allocator));
    o2 = double_to_sign(predicates::orient2d(
        &v1[1..],
        &v2[1..],
        &v3[1..],
        allocator,
    ));
    let oo4 = o2;
    if o1 != o2 {
        return false;
    };

    // Projection on (x,z)-plane -> p VS v1
    let pxz = [p[0], p[2]];
    let v1xz = [v1[0], v1[2]];
    let v2xz = [v2[0], v2[2]];
    let v3xz = [v3[0], v3[2]];
    o1 = double_to_sign(predicates::orient2d(&pxz, &v2xz, &v3xz, allocator));
    o2 = double_to_sign(predicates::orient2d(&v1xz, &v2xz, &v3xz, allocator));
    let oo6 = o2;
    if o1 != o2 {
        return false;
    };

    // Projection on (x,y)-plane -> p VS v2
    o1 = double_to_sign(predicates::orient2d(p, v3, v1, allocator));
    o2 = oo2;
    if o1 != o2 {
        return false;
    };

    // Projection on (y,z)-plane -> p VS v2
    o1 = double_to_sign(predicates::orient2d(&p[1..], &v3[1..], &v1[1..], allocator));
    o2 = oo4;
    if o1 != o2 {
        return false;
    };

    // Projection on (x,z)-plane -> p VS v2
    o1 = double_to_sign(predicates::orient2d(&pxz, &v3xz, &v1xz, allocator));
    o2 = oo6;
    if o1 != o2 {
        return false;
    };

    // Projection on (x,y)-plane -> p VS v3
    o1 = double_to_sign(predicates::orient2d(p, v1, v2, allocator));
    o2 = oo2;
    if o1 != o2 {
        return false;
    };

    // Projection on (y,z)-plane -> p VS v3
    o1 = double_to_sign(predicates::orient2d(&p[1..], &v1[1..], &v2[1..], allocator));
    o2 = oo4;
    if o1 != o2 {
        return false;
    };

    // Projection on (x,z)-plane -> p VS v3
    o1 = double_to_sign(predicates::orient2d(&pxz, &v1xz, &v2xz, allocator));
    o2 = oo6;
    if o1 != o2 {
        return false;
    };

    true
}

#[inline]
pub fn inner_segment_cross_inner_triangle<A: Allocator + Copy>(
    u1: &[f64],
    u2: &[f64],
    v1: &[f64],
    v2: &[f64],
    v3: &[f64],
    allocator: A,
) -> bool {
    let mut bound = u1[0].min(u2[0]); // min(u1,u2) alogng x-axis
    if v1[0] <= bound && v2[0] <= bound && v3[0] <= bound {
        return false;
    }
    bound = u1[0].max(u2[0]); // max(u1,u2) alogng x-axis
    if v1[0] >= bound && v2[0] >= bound && v3[0] >= bound {
        return false;
    }
    bound = u1[1].min(u2[1]); // min(u1,u2) alogng y-axis
    if v1[1] <= bound && v2[1] <= bound && v3[1] <= bound {
        return false;
    }
    bound = u1[1].max(u2[1]); // max(u1,u2) alogng y-axis
    if v1[1] >= bound && v2[1] >= bound && v3[1] >= bound {
        return false;
    }
    bound = u1[2].min(u2[2]); // min(u1,u2) alogng z-axis
    if v1[2] <= bound && v2[2] <= bound && v3[2] <= bound {
        return false;
    }
    bound = u1[2].max(u2[2]); // max(u1,u2) alogng z-axis
    if v1[2] >= bound && v2[2] >= bound && v3[2] >= bound {
        return false;
    }

    let orient_u1_tri = double_to_sign(predicates::orient3d(u1, v1, v2, v3, allocator));
    let orient_u2_tri = double_to_sign(predicates::orient3d(u2, v1, v2, v3, allocator));

    // Check if triangle vertices and at least one of the segment endpoints are coplanar:
    // in this case there is no proper intersection.
    if orient_u1_tri == Orientation::Zero || orient_u2_tri == Orientation::Zero {
        return false;
    }

    // Endpoints of one segment cannot stay both in one of the same half-space defined by the triangle.
    if orient_u1_tri == orient_u2_tri {
        return false;
    }

    // Since now, endpoints are one abouve and one below the triangle-plane.

    // Intersection between segment and triangle sides are not proper.
    // Check also if segment intersect the triangle-plane outside the triangle.
    let orient_u_v1v2 = double_to_sign(predicates::orient3d(u1, u2, v1, v2, allocator));
    let orient_u_v2v3 = double_to_sign(predicates::orient3d(u1, u2, v2, v3, allocator));

    if orient_u_v1v2 == Orientation::Zero || orient_u_v2v3 == Orientation::Zero {
        return false;
    }

    if orient_u_v1v2 != orient_u_v2v3 {
        return false;
    }

    let orient_u_v3v1 = double_to_sign(predicates::orient3d(u1, u2, v3, v1, allocator));

    if orient_u_v3v1 == Orientation::Zero {
        return false;
    }
    if orient_u_v3v1 != orient_u_v2v3 {
        return false;
    }

    // Finally, we have a proper intersection.
    return true;
}

#[inline(always)]
pub fn same_point(p: &[f64], q: &[f64]) -> bool {
    p[0] == q[0] && p[1] == q[1] && p[2] == q[2]
}

pub fn same_half_plane<A: Allocator + Copy>(
    p: &[f64],
    q: &[f64],
    v1: &[f64],
    v2: &[f64],
    allocator: A,
) -> bool {
    // Projection on (x,y)-plane
    if double_to_sign(predicates::orient2d(p, v1, v2, allocator))
        != double_to_sign(predicates::orient2d(q, v1, v2, allocator))
    {
        return false;
    }

    // Projection on (y,z)-plane
    if double_to_sign(predicates::orient2d(&p[1..], &v1[1..], &v2[1..], allocator))
        != double_to_sign(predicates::orient2d(&q[1..], &v1[1..], &v2[1..], allocator))
    {
        return false;
    }

    // Projection on (x,z)-plane
    let pxz = [p[0], p[2]];
    let qxz = [q[0], q[2]];
    let v1xz = [v1[0], v1[2]];
    let v2xz = [v2[0], v2[2]];
    return double_to_sign(predicates::orient2d(&pxz, &v1xz, &v2xz, allocator))
        == double_to_sign(predicates::orient2d(&qxz, &v1xz, &v2xz, allocator));
}

#[inline]
pub fn mis_alignment<A: Allocator + Copy>(p: &[f64], q: &[f64], r: &[f64], allocator: A) -> bool {
    // Projection on (x,y)-plane
    if predicates::orient2d(p, q, r, allocator) != 0.0 {
        return true;
    }

    // Projection on (y,z)-plane
    if predicates::orient2d(&p[1..], &q[1..], &r[1..], allocator) != 0.0 {
        return true;
    }

    // Projection on (x,z)-plane
    let pxz = [p[0], p[2]];
    let qxz = [q[0], q[2]];
    let rxz = [r[0], r[2]];
    return predicates::orient2d(&pxz, &qxz, &rxz, allocator) != 0.0;
}

#[inline]
pub fn inner_segments_cross<A: Allocator + Copy>(
    u1: &[f64],
    u2: &[f64],
    v1: &[f64],
    v2: &[f64],
    allocator: A,
) -> bool {
    // The 4 endpoints must be coplanar
    if predicates::orient3d(u1, u2, v1, v2, allocator) != 0.0 {
        return false;
    }

    // Endpoints of one segment cannot stay either on the same side of the other one.
    if same_half_plane(u1, u2, v1, v2, allocator) || same_half_plane(v1, v2, u1, u2, allocator) {
        return false;
    }

    // Each segment endpoint cannot be aligned with the other segment.
    if !mis_alignment(u1, v1, v2, allocator) {
        return false;
    };
    if !mis_alignment(u2, v1, v2, allocator) {
        return false;
    };
    if !mis_alignment(v1, u1, u2, allocator) {
        return false;
    };
    if !mis_alignment(v2, u1, u2, allocator) {
        return false;
    };

    // If the segment projected on one coordinate plane cross -> segmant cross.
    // Projection on (x,y)-plane
    if predicates::orient2d(u1, u2, v1, allocator) != 0.0 {
        return true;
    };
    if predicates::orient2d(v1, v2, u2, allocator) != 0.0 {
        return true;
    };
    // Projection on (y,z)-plane
    if predicates::orient2d(&u1[1..], &u2[1..], &v1[1..], allocator) != 0.0 {
        return true;
    };
    if predicates::orient2d(&v1[1..], &v2[1..], &u2[1..], allocator) != 0.0 {
        return true;
    };
    // Projection on (z,x)-plane
    let u1xz = [u1[0], u1[2]];
    let u2xz = [u2[0], u2[2]];
    let v1xz = [v1[0], v1[2]];
    let v2xz = [v2[0], v2[2]];
    if predicates::orient2d(&u1xz, &u2xz, &v1xz, allocator) != 0.0 {
        return true;
    };
    if predicates::orient2d(&v1xz, &v2xz, &u2xz, allocator) != 0.0 {
        return true;
    };

    return false;
}

#[inline(always)]
pub fn point_in_inner_segment<A: Allocator + Copy>(
    p: &[f64],
    v1: &[f64],
    v2: &[f64],
    allocator: A,
) -> bool {
    return !mis_alignment(p, v1, v2, allocator)
        && ((v1[0] < v2[0] && v1[0] < p[0] && p[0] < v2[0])
            || (v1[0] > v2[0] && v1[0] > p[0] && p[0] > v2[0])
            || (v1[1] < v2[1] && v1[1] < p[1] && p[1] < v2[1])
            || (v1[1] > v2[1] && v1[1] > p[1] && p[1] > v2[1])
            || (v1[2] < v2[2] && v1[2] < p[2] && p[2] < v2[2])
            || (v1[2] > v2[2] && v1[2] > p[2] && p[2] > v2[2]));
}

#[inline(always)]
pub fn point_in_segment<A: Allocator + Copy>(
    p: &[f64],
    v1: &[f64],
    v2: &[f64],
    allocator: A,
) -> bool {
    return same_point(p, v1) || same_point(p, v2) || point_in_inner_segment(p, v1, v2, allocator);
}

/// is the point p in the segment v1v2?
/// assuming that p, v1, v2 are collinear
#[inline(always)]
pub fn point_in_segment_general<A: Allocator + Copy>(
    p: &Point3D,
    v1: &Point3D,
    v2: &Point3D,
    allocator: A
) -> bool {
    let comp1 = less_than_on_x(v1, p, allocator);
    let comp2 = less_than_on_x(p, v2, allocator);
    if comp1 != comp2 {
        if comp1 == Orientation::Zero || comp2 == Orientation::Zero {
            return true;
        } else {
            return false;
        }
    } else if comp1 != Orientation::Zero {
        return true;
    }

    let comp1 = less_than_on_y(v1, p, allocator);
    let comp2 = less_than_on_y(p, v2, allocator);

    if comp1 != comp2 {
        if comp1 == Orientation::Zero || comp2 == Orientation::Zero {
            return true;
        } else {
            return false;
        }
    } else if comp1 != Orientation::Zero {
        return true;
    }

    let comp1 = less_than_on_z(v1, p, allocator);
    let comp2 = less_than_on_z(p, v2, allocator);
    if comp1 != comp2 {
        if comp1 == Orientation::Zero || comp2 == Orientation::Zero {
            return true;
        } else {
            return false;
        }
    } else {
        // the segment is a point
        return true;
    }
}

#[inline(always)]
pub fn inner_segment_cross_triangle<A: Allocator + Copy>(
    u1: &[f64],
    u2: &[f64],
    v1: &[f64],
    v2: &[f64],
    v3: &[f64],
    allocator: A,
) -> bool {
    return point_in_inner_segment(&v1, &u1, &u2, allocator)
        || point_in_inner_segment(&v2, &u1, &u2, allocator)
        || point_in_inner_segment(&v3, &u1, &u2, allocator)
        || inner_segments_cross(&v2, &v3, &u1, &u2, allocator)
        || inner_segments_cross(&v3, &v1, &u1, &u2, allocator)
        || inner_segments_cross(&v1, &v2, &u1, &u2, allocator)
        || inner_segment_cross_inner_triangle(&u1, &u2, &v1, &v2, &v3, allocator);
}
#[inline(always)]
pub fn point_in_triangle<A: Allocator + Copy>(
    p: &[f64],
    v1: &[f64],
    v2: &[f64],
    v3: &[f64],
    allocator: A,
) -> bool {
    return point_in_segment(p, v1, v2, allocator)
        || point_in_segment(p, v2, v3, allocator)
        || point_in_segment(p, v3, v1, allocator)
        || point_in_inner_triangle(p, v1, v2, v3, allocator);
}

#[inline(always)]
pub fn orient2d_by_axis<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    axis: usize,
    bump: A,
) -> Orientation {
    let ori = if axis == 0 {
        orient2d(&pa[1..], &pb[1..], &pc[1..], bump)
    } else if axis == 1 {
        let pa = [pa[2], pa[0]];
        let pb = [pb[2], pb[0]];
        let pc = [pc[2], pc[0]];
        orient2d(&pa, &pb, &pc, bump)
    } else {
        orient2d(pa, pb, pc, bump)
    };
    double_to_sign(ori)
}
