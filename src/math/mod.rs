use std::ops::{Add, Mul, Sub};

#[inline(always)]
pub fn sub_in<T: Copy + Sub<Output = T>>(a: &[T], b: &[T], c: &mut [T]) {
    for i in 0..a.len() {
        c[i] = a[i] - b[i];
    }
}

#[inline(always)]
pub fn sub_short<const N: usize, T: Copy + Sub<Output = T>>(a: &[T], b: &[T]) -> [T; N] {
    let mut c = [a[0] - b[0]; N];
    for i in 1..N {
        c[i] = a[i] - b[i];
    }
    c
}

#[inline(always)]
pub fn norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |acc, x| acc + x * x).sqrt()
}

#[inline(always)]
pub fn square_norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |acc, x| acc + x * x)
}

#[inline(always)]
pub fn dot<T: Copy + Add<Output = T> + Mul<Output = T>>(a: &[T], b: &[T]) -> T {
    let mut ret = a[0] * b[0];
    for i in 1..a.len() {
        ret = ret + a[i] * b[i];
    }
    ret
}

#[inline(always)]
pub fn cross_in<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>>(
    a: &[T],
    b: &[T],
    c: &mut [T],
) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

#[inline(always)]
pub fn cross<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>>(
    a: &[T],
    b: &[T],
) -> [T; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
