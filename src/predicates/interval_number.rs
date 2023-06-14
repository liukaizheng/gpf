use core::f64;
use std::ops::{Add, Mul, Sub};

#[derive(Clone)]
pub struct IntervalNumber {
    low: f64,
    high: f64,
}

impl IntervalNumber {
    #[inline(always)]
    pub fn new(inf: f64, sup: f64) -> Self {
        Self {
            low: inf,
            high: sup,
        }
    }
}

impl Default for IntervalNumber {
    #[inline(always)]
    fn default() -> Self {
        Self {
            low: f64::NAN,
            high: f64::NAN,
        }
    }
}

impl From<f64> for IntervalNumber {
    #[inline(always)]
    fn from(val: f64) -> Self {
        Self {
            low: -val,
            high: val,
        }
    }
}

impl IntervalNumber {
    #[inline(always)]
    pub fn not_zero(&self) -> bool {
        self.low < 0.0 || self.high < 0.0
    }

    #[inline(always)]
    pub fn positive(&self) -> bool {
        self.low < 0.0
    }

    #[inline(always)]
    pub fn negative(&self) -> bool {
        self.high < 0.0
    }

    #[inline(always)]
    pub fn neg(&mut self) {
        std::mem::swap(&mut self.low, &mut self.high);
    }
}

impl PartialEq for IntervalNumber {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.low == other.low && self.high == other.high
    }
}

#[inline(always)]
fn add(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    IntervalNumber {
        low: (a.low + b.low).next_up(),
        high: (a.high + b.high).next_up(),
    }
}

#[inline(always)]
fn sub(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    IntervalNumber {
        low: (a.low + b.high).next_up(),
        high: (a.high + b.low).next_up(),
    }
}

fn mul(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    let signs: (bool, bool, bool, bool) = (
        a.low.is_sign_negative() & (a.low != 0.0),
        a.high.is_sign_negative() & (a.high != 0.0),
        b.low.is_sign_negative() & (b.low != 0.0),
        b.high.is_sign_negative() & (b.high != 0.0),
    );
    match signs {
        (true, false, true, false) => {
            IntervalNumber::new((a.low * (-b.low)).next_up(), (a.high * b.high).next_up())
        }
        (true, false, false, true) => {
            IntervalNumber::new((a.high * b.low).next_up(), ((-a.low) * b.high).next_up())
        }
        (true, false, false, false) => {
            IntervalNumber::new((a.high * b.low).next_up(), (a.high * b.high).next_up())
        }
        (false, true, true, false) => {
            IntervalNumber::new((a.low * b.high).next_up(), (a.high * (-b.low)).next_up())
        }
        (false, true, false, true) => {
            IntervalNumber::new(((-a.high) * b.high).next_up(), (a.low * b.low).next_up())
        }
        (false, true, false, false) => {
            IntervalNumber::new((a.low * b.high).next_up(), (a.low * b.low).next_up())
        }
        (false, false, true, false) => {
            IntervalNumber::new((a.low * b.high).next_up(), (a.high * b.high).next_up())
        }
        (false, false, false, true) => {
            IntervalNumber::new((a.high * b.low).next_up(), (a.low * b.low).next_up())
        }
        (false, false, false, false) => IntervalNumber::new(
            ((a.low * b.high).max(a.high * b.low)).next_up(),
            ((a.low * b.low).max(a.high * b.high)).next_up(),
        ),
        _ => IntervalNumber::default(),
    }
}

macro_rules! impl_op {
    (trait $op: ident, $func: ident) => {
        impl $op for IntervalNumber {
            type Output = Self;
            #[inline(always)]
            fn $func(self, other: Self) -> Self {
                $func(&self, &other)
            }
        }

        impl $op<&IntervalNumber> for IntervalNumber {
            type Output = Self;
            #[inline(always)]
            fn $func(self, other: &Self) -> Self {
                $func(&self, other)
            }
        }

        impl $op<IntervalNumber> for &IntervalNumber {
            type Output = IntervalNumber;
            #[inline(always)]
            fn $func(self, other: IntervalNumber) -> IntervalNumber {
                $func(self, &other)
            }
        }

        impl $op<&IntervalNumber> for &IntervalNumber {
            type Output = IntervalNumber;
            #[inline(always)]
            fn $func(self, other: &IntervalNumber) -> IntervalNumber {
                $func(self, other)
            }
        }
    };
}

impl_op!(trait Add, add);
impl_op!(trait Sub, sub);
impl_op!(trait Mul, mul);

#[test]
fn interval_number_operation() {
    let a = IntervalNumber::from(2.0);
    let b = &a - &a - &a;
    assert!(b.high.round() == -2.0);
}
