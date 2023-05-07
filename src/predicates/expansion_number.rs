use std::ops::{Add, Deref, Sub, Mul};

use bumpalo::{collections::vec::Vec, Bump};

use super::{fast_expansion_sum_zeroelim, fast_expansion_diff_zeroelim, mul_expansion_zeroelim};

#[derive(Clone)]
struct ExpansionNumOwn<'b> {
    vec: Vec<'b, f64>,
}

#[derive(Clone)]
struct ExpansionNumRef<'a, 'b> {
    vec: &'a Vec<'b, f64>,
}

impl<'b> Deref for ExpansionNumOwn<'b> {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.vec
    }
}

impl<'a, 'b> Deref for ExpansionNumRef<'a, 'b> {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.vec
    }
}

pub trait ExpansionNum<'b>: Deref<Target = [f64]> {
    fn bump(&self) -> &'b Bump;
}

impl<'b> ExpansionNum<'b> for ExpansionNumOwn<'b> {
    fn bump(&self) -> &'b Bump {
        self.vec.bump()
    }
}

impl<'a, 'b> ExpansionNum<'b> for ExpansionNumRef<'a, 'b> {
    fn bump(&self) -> &'b Bump {
        self.vec.bump()
    }
}

#[inline]
fn add<'b, T1: ExpansionNum<'b>, T2: ExpansionNum<'b>>(t1: T1, t2: T2) -> Vec<'b, f64> {
    fast_expansion_sum_zeroelim(&t1, &t2, t1.bump())
}

#[inline]
fn sub<'b, T1: ExpansionNum<'b>, T2: ExpansionNum<'b>>(t1: T1, t2: T2) -> Vec<'b, f64> {
    fast_expansion_diff_zeroelim(&t1, &t2, t1.bump())
}

#[inline]
fn mul<'b, T1: ExpansionNum<'b>, T2: ExpansionNum<'b>>(t1: T1, t2: T2) -> Vec<'b, f64> {
    mul_expansion_zeroelim(&t1, &t2, t1.bump())
}

macro_rules! impl_op {
    (trait $op: ident, $func: ident) => {
        impl<'b, RHS: ExpansionNum<'b>> $op<RHS> for ExpansionNumOwn<'b> {
            type Output = ExpansionNumOwn<'b>;

            #[inline]
            fn $func(self, rhs: RHS) -> Self::Output {
                Self::Output {
                    vec: $func(self, rhs),
                }
            }
        }

        impl<'b, RHS: ExpansionNum<'b>> $op<RHS> for ExpansionNumRef<'_, 'b> {
            type Output = ExpansionNumOwn<'b>;

            #[inline]
            fn $func(self, rhs: RHS) -> Self::Output {
                Self::Output {
                    vec: $func(self, rhs),
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
    let v1 = bumpalo::vec![in &bump; 2.0; 1];
    let v2 = bumpalo::vec![in &bump; 1.0; 1];
    let v1_own = ExpansionNumOwn { vec: v1 };
    {
        let v2_ref = ExpansionNumRef { vec: &v2 };
        let v3 = v1_own.clone() + v2_ref;
        assert_eq!(v3[0], 3.0);
    }
    let v4 = v1_own + ExpansionNumOwn { vec: v2 };
    assert_eq!(v4[0], 3.0);
}
