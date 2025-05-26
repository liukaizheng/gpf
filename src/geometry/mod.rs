mod bbox;
mod curve;
mod surface;

pub use bbox::BBox;
pub use curve::*;
pub use surface::*;

use crate::math::{cross, dot, square_norm};

#[inline]
fn is_parallel(a: &[f64], b: &[f64], cos_tol: f64) -> bool {
    return (dot(a, b).abs() - 1.0) < cos_tol;
}

#[inline]
fn sq_dist_to_line(p: &[f64], o: &[f64], d: &[f64]) -> f64 {
    let v = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
    square_norm(&cross(&v, d))
}