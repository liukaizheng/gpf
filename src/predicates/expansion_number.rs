use std::ops::{Add, Deref, Sub, Mul};

use bumpalo::{collections::vec::Vec, Bump};

use super::{predicates::{fast_expansion_sum_zeroelim, fast_expansion_diff_zeroelim, mul_expansion_zeroelim}};

#[derive(Clone)]
pub struct ExpansionNum<'b> {
    vec: Vec<'b, f64>,
}

impl <'b> ExpansionNum<'b> {
    #[inline]
    pub fn bump(&self) -> &'b Bump {
        self.vec.bump()
    }
}


impl<'b> Deref for ExpansionNum<'b> {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.vec
    }
}

impl <'b> From<Vec<'b, f64>> for ExpansionNum<'b> {
    #[inline(always)]
    fn from(vec: Vec<'b, f64>) -> Self {
        Self {
            vec
        }
    }
}


#[inline]
fn add<'b>(a: &[f64], b: &[f64], bump: &'b Bump) -> Vec<'b, f64> {
    fast_expansion_sum_zeroelim(a, b, bump)
}

#[inline]
fn sub<'b>(a: &[f64], b: &[f64], bump: &'b Bump) -> Vec<'b, f64> {
    fast_expansion_diff_zeroelim(a, b, bump)
}
#[inline]
fn mul<'b>(a: &[f64], b: &[f64], bump: &'b Bump) -> Vec<'b, f64> {
    mul_expansion_zeroelim(a, b, bump)
}

macro_rules! impl_op {
    (trait $op: ident, $func: ident) => {
        impl<'b> $op for ExpansionNum<'b> {
            type Output = ExpansionNum<'b>;

            #[inline]
            fn $func(self, rhs: Self) -> Self::Output {
                Self::Output {
                    vec: $func(&self, &rhs, self.bump()),
                }
            }
        }

        impl<'b> $op for &ExpansionNum<'b> {
            type Output = ExpansionNum<'b>;

            #[inline]
            fn $func(self, rhs: Self) -> Self::Output {
                Self::Output {
                    vec: $func(self, rhs, self.bump()),
                }
            }
        }

        impl<'b> $op<&ExpansionNum<'b>> for ExpansionNum<'b> {
            type Output = ExpansionNum<'b>;

            #[inline]
            fn $func(self, rhs: &Self) -> Self::Output {
                Self::Output {
                    vec: $func(&self, rhs, self.bump()),
                }
            }
        }

        impl<'b> $op<ExpansionNum<'b>> for &ExpansionNum<'b> {
            type Output = ExpansionNum<'b>;

            #[inline]
            fn $func(self, rhs: ExpansionNum<'b>) -> Self::Output {
                Self::Output {
                    vec: $func(self, &rhs, self.bump()),
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
    let bump = Bump::new();
    let v1 = ExpansionNum{ vec: bumpalo::vec![in &bump; 2.0; 1]};
    let v2 = ExpansionNum { vec: bumpalo::vec![in &bump; 1.0; 1] };
    let v3 = v1 + &v2;
    assert_eq!(v3[0], 3.0);
    let v4 = v3 + v2;
    assert_eq!(v4[0], 4.0);
}
