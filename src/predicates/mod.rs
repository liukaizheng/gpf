mod expansion_number;
mod generic_point;
mod interval_number;
mod orient2d;
mod predicates;

use std::ops::{Add, Mul, Sub};

pub use expansion_number::*;
pub use generic_point::*;
pub use interval_number::*;
pub use orient2d::*;

pub enum Orientation {
    Positive,
    Negative,
    Zero,
}

#[inline(always)]
fn double_to_sign(x: f64) -> Orientation {
    if x > 0.0 {
        Orientation::Positive
    } else if x < 0.0 {
        Orientation::Negative
    } else {
        Orientation::Zero
    }
}
