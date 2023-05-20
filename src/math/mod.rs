use std::ops::{Add, Mul, Sub};

#[inline(always)]
pub fn sub<T: Copy + Sub<Output = T>>(a: &[T], b: &[T], c: &mut [T]) {
    for i in 0..a.len() {
        c[i] = a[i] - b[i];
    }
}

#[inline(always)]
pub fn norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |acc, x| acc + x * x).sqrt()
}

#[inline(always)]
pub fn cross<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>>(
    a: &[T],
    b: &[T],
    c: &mut [T],
) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}