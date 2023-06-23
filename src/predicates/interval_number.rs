use core::f64;
use std::ops::{Add, Mul, Sub};
use std::simd::{f64x2, u64x2, SimdFloat, ToBitMask};

const SIGN_LOW_MASK: u64x2 = u64x2::from_array([0_u64, (1_u64 << 63)]);
const SIGN_HIGH_MASK: u64x2 = u64x2::from_array([(1_u64 << 63), 0u64]);

#[derive(Clone)]
pub struct IntervalNumber {
    /// [high; low]
    data: f64x2,
}

impl IntervalNumber {
    #[inline(always)]
    pub fn new(low: f64, high: f64) -> Self {
        Self {
            data: f64x2::from_array([high, low]),
        }
    }
}

impl Default for IntervalNumber {
    #[inline(always)]
    fn default() -> Self {
        Self {
            data: f64x2::from_array([f64::NAN, f64::NAN]),
        }
    }
}

impl From<f64> for IntervalNumber {
    #[inline(always)]
    fn from(val: f64) -> Self {
        Self {
            data: f64x2::from_array([val, -val]),
        }
    }
}

impl IntervalNumber {
    #[inline(always)]
    pub fn not_zero(&self) -> bool {
        self.data[0] < 0.0 || self.data[1] < 0.0
    }

    #[inline(always)]
    pub fn positive(&self) -> bool {
        self.data[0] < 0.0
    }

    #[inline(always)]
    pub fn negative(&self) -> bool {
        self.data[1] < 0.0
    }

    #[inline(always)]
    pub fn neg(&mut self) {
        self.data = self.data.reverse();
    }

    pub fn from_f64x2(data: f64x2) -> Self {
        Self { data }
    }
}

impl PartialEq for IntervalNumber {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

#[inline(always)]
fn next_up(data: &mut f64x2) {
    data[0] = data[0].next_up();
    data[1] = data[1].next_up();
}

#[inline(always)]
fn add(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    let mut data = a.data + b.data;
    next_up(&mut data);
    IntervalNumber { data }
}

#[inline(always)]
fn sub(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    let mut data = b.data.reverse() + a.data;
    next_up(&mut data);
    IntervalNumber { data }
}

fn mul(a: &IntervalNumber, b: &IntervalNumber) -> IntervalNumber {
    let mask =
        (a.data.is_sign_negative().to_bitmask() << 2) | (b.data.is_sign_negative().to_bitmask());
    let mut data = match mask {
        0 => {
            let llhh = a.data * b.data;
            let lhhl = a.data * b.data.reverse();
            [lhhl.reduce_max(), llhh.reduce_max()].into()
        }

        1 => a.data.reverse() * f64x2::from([b.data[1], b.data[1]]),

        2 => a.data * f64x2::from([b.data[0], b.data[0]]),

        4 => f64x2::from([a.data[1], a.data[1]]) * b.data.reverse(),

        5 => {
            let mut ai = a.data.to_bits();
            ai ^= SIGN_HIGH_MASK;
            (f64x2::from_bits(ai) * b.data).reverse()
        }

        6 => {
            let mut bi = b.data.to_bits();
             bi ^= SIGN_LOW_MASK;
            a.data * f64x2::from_bits(bi).reverse()
        }

        8 => (f64x2::from([a.data[0]; 2]) * b.data).reverse(),

        9 => {
            let mut ai = a.data.to_bits();
            ai ^= SIGN_LOW_MASK;
            f64x2::from_bits(ai).reverse() * b.data
        }
        10 => {
            let mut bi = b.data.to_bits();
            bi ^= SIGN_LOW_MASK;
            a.data * f64x2::from_bits(bi)
        }

        _ => {
            panic!("never")
        }
    };
    next_up(&mut data);
    IntervalNumber { data }
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
    assert!(b.data[0].round() == -2.0);
}
