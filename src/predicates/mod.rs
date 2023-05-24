mod expansion_number;
mod generic_point;
mod interval_number;
mod orient2d;
pub mod orient3d;
mod predicates;

use bumpalo::{collections::Vec, Bump};
use std::ops::{Add, Mul, Sub};

pub use expansion_number::*;
pub use generic_point::*;
pub use interval_number::*;
pub use orient2d::*;
pub use predicates::*;

#[derive(PartialEq, Eq)]
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
fn abs_max(vec: Vec<'_, f64>) -> Option<f64> {
    vec.into_iter().map(|e| e.abs()).reduce(|acc, e| acc.max(e))
}

#[inline(always)]
fn dummy_abs_max<T>(_: Vec<'_, T>) -> Option<T> {
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
pub fn get_exponent(x: f64) -> i32 {
    if x == 0.0 {
        0
    } else {
        x.to_bits().wrapping_shr(52) as i32 - 1023
    }
}

fn max_comp_in_tri_normal_impl<
    'b,
    const NEED_MAX: bool,
    T: 'b + GenericNum,
    F: FnOnce(Vec<'b, T>) -> Option<T>,
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
    bump: &'b Bump,
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
        abs_max(bumpalo::vec![in bump;v3x, v3y, v3z, v2x, v2y, v2z])
    } else {
        None
    };
    (nvx, nvy, nvz, max_var)
}

#[inline(always)]
fn max_comp_in_tri_normal_filter<'b>(
    ov1x: f64,
    ov1y: f64,
    ov1z: f64,
    ov2x: f64,
    ov2y: f64,
    ov2z: f64,
    ov3x: f64,
    ov3y: f64,
    ov3z: f64,
    bump: &'b Bump,
) -> usize {
    let (nvx, nvy, nvz, max_var) = max_comp_in_tri_normal_impl::<'b, true, _, _>(
        ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z, abs_max, bump,
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
fn max_comp_in_tri_normal_exact<'b>(
    ov1x: ExpansionNum<'b>,
    ov1y: ExpansionNum<'b>,
    ov1z: ExpansionNum<'b>,
    ov2x: ExpansionNum<'b>,
    ov2y: ExpansionNum<'b>,
    ov2z: ExpansionNum<'b>,
    ov3x: ExpansionNum<'b>,
    ov3y: ExpansionNum<'b>,
    ov3z: ExpansionNum<'b>,
) -> usize {
    let bump = ov1x.bump();
    let (nvx, nvy, nvz, _) = max_comp_in_tri_normal_impl::<'b, false, _, _>(
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
        bump,
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
pub fn max_comp_in_tri_normal<'b>(ov1: &[f64], ov2: &[f64], ov3: &[f64], bump: &'b Bump) -> usize {
    let max_comp = max_comp_in_tri_normal_filter(
        ov1[0], ov1[1], ov1[2], ov2[0], ov2[1], ov2[2], ov3[0], ov3[1], ov3[2], bump,
    );
    if max_comp < 4 {
        max_comp
    } else {
        max_comp_in_tri_normal_exact(
            bumpalo::vec![in bump;ov1[0]].into(),
            bumpalo::vec![in bump;ov1[1]].into(),
            bumpalo::vec![in bump;ov1[2]].into(),
            bumpalo::vec![in bump;ov2[0]].into(),
            bumpalo::vec![in bump;ov2[1]].into(),
            bumpalo::vec![in bump;ov2[2]].into(),
            bumpalo::vec![in bump;ov3[0]].into(),
            bumpalo::vec![in bump;ov3[1]].into(),
            bumpalo::vec![in bump;ov3[2]].into(),
        )
    }
}
