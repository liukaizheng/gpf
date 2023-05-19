mod expansion_number;
mod generic_point;
mod interval_number;
mod orient2d;
pub mod orient3d;
mod predicates;

use bumpalo::collections::Vec;
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
