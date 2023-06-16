use std::{
    alloc::{Allocator, Global},
    ops::{Add, Deref, DerefMut, Mul, Sub},
};

use super::predicates::{
    fast_expansion_diff_zeroelim, fast_expansion_sum_zeroelim, mul_expansion_zeroelim,
};

#[derive(Clone)]
pub struct ExpansionNum<A: Allocator + Copy = Global> {
    pub vec: Vec<f64, A>,
}

impl<A: Allocator + Copy> ExpansionNum<A> {
    #[inline(always)]
    pub fn allocator(&self) -> &A {
        self.vec.allocator()
    }
}

impl<A: Allocator + Copy> Deref for ExpansionNum<A> {
    type Target = [f64];

    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        &self.vec
    }
}

impl<A: Allocator + Copy> DerefMut for ExpansionNum<A> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.vec
    }
}

impl<A: Allocator + Copy> From<Vec<f64, A>> for ExpansionNum<A> {
    #[inline(always)]
    fn from(vec: Vec<f64, A>) -> Self {
        Self { vec }
    }
}

impl<A: Allocator + Copy> ExpansionNum<A> {
    #[inline(always)]
    pub fn not_zero(&self) -> bool {
        *self.last().unwrap() != 0.0
    }
    #[inline(always)]
    pub fn negative(&self) -> bool {
        *self.last().unwrap() < 0.0
    }

    #[inline(always)]
    pub fn neg(&mut self) {
        for v in &mut self.vec {
            *v = -*v;
        }
    }
}

#[inline(always)]
fn add<A: Allocator + Copy>(a: &[f64], b: &[f64], allocator: A) -> Vec<f64, A> {
    fast_expansion_sum_zeroelim(a, b, allocator)
}

#[inline(always)]
fn sub<A: Allocator + Copy>(a: &[f64], b: &[f64], allocator: A) -> Vec<f64, A> {
    fast_expansion_diff_zeroelim(a, b, allocator)
}

#[inline(always)]
fn mul<A: Allocator + Copy>(a: &[f64], b: &[f64], allocator: A) -> Vec<f64, A> {
    mul_expansion_zeroelim(a, b, allocator)
}

macro_rules! impl_op {
    (trait $op: ident, $func: ident) => {
        impl<A: Allocator + Copy> $op for ExpansionNum<A> {
            type Output = ExpansionNum<A>;

            #[inline(always)]
            fn $func(self, rhs: Self) -> Self::Output {
                Self::Output {
                    vec: $func(&self, &rhs, *self.allocator()),
                }
            }
        }

        impl<A: Allocator + Copy> $op for &ExpansionNum<A> {
            type Output = ExpansionNum<A>;

            #[inline(always)]
            fn $func(self, rhs: Self) -> Self::Output {
                Self::Output {
                    vec: $func(self, rhs, *self.allocator()),
                }
            }
        }

        impl<A: Allocator + Copy> $op<&ExpansionNum<A>> for ExpansionNum<A> {
            type Output = ExpansionNum<A>;

            #[inline(always)]
            fn $func(self, rhs: &Self) -> Self::Output {
                Self::Output {
                    vec: $func(&self, rhs, *self.allocator()),
                }
            }
        }

        impl<A: Allocator + Copy> $op<ExpansionNum<A>> for &ExpansionNum<A> {
            type Output = ExpansionNum<A>;

            #[inline(always)]
            fn $func(self, rhs: ExpansionNum<A>) -> Self::Output {
                Self::Output {
                    vec: $func(self, &rhs, *self.allocator()),
                }
            }
        }
    };
}

impl_op!(trait Add, add);
impl_op!(trait Sub, sub);
impl_op!(trait Mul, mul);

#[test]
fn test_expansion_operations() {
    let v1 = ExpansionNum { vec: vec![2.0; 1] };
    let v2 = ExpansionNum { vec: vec![1.0; 1] };
    let v3 = v1 + &v2;
    assert_eq!(v3[0], 3.0);
    let v4 = v3 + v2;
    assert_eq!(v4[0], 4.0);
}
