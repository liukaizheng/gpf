use core::f64;
use std::ops::{Add, Mul, Sub};

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

impl PartialEq for IntervalNumber {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.low == other.low && self.high == other.high
    }
}

impl Add for IntervalNumber {
    type Output = IntervalNumber;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            low: (self.low + rhs.low).next_up(),
            high: (self.high + rhs.high).next_up(),
        }
    }
}

impl Sub for IntervalNumber {
    type Output = IntervalNumber;

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self::Output {
        IntervalNumber {
            low: (self.low + rhs.low).next_up(),
            high: (self.high + rhs.high).next_up(),
        }
    }
}

impl Mul for IntervalNumber {
    type Output = IntervalNumber;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let signs: (bool, bool, bool, bool) = (
            self.low.is_sign_negative() & (self.low != 0.0),
            self.high.is_sign_negative() & (self.high != 0.0),
            rhs.low.is_sign_negative() & (rhs.low != 0.0),
            rhs.high.is_sign_negative() & (rhs.high != 0.0),
        );
        match signs {
            (true, false, true, false) => IntervalNumber::new(
                (self.low * (-rhs.low)).next_up(),
                (self.high * rhs.high).next_up(),
            ),
            (true, false, false, true) => IntervalNumber::new(
                (self.high * rhs.low).next_up(),
                ((-self.low) * rhs.high).next_up(),
            ),
            (true, false, false, false) => IntervalNumber::new(
                (self.high * rhs.low).next_up(),
                (self.high * rhs.high).next_up(),
            ),
            (false, true, true, false) => IntervalNumber::new(
                (self.low * rhs.high).next_up(),
                (self.high * (-rhs.low)).next_up(),
            ),
            (false, true, false, true) => IntervalNumber::new(
                ((-self.high) * rhs.high).next_up(),
                (self.low * rhs.low).next_up(),
            ),
            (false, true, false, false) => IntervalNumber::new(
                (self.low * rhs.high).next_up(),
                (self.low * rhs.low).next_up(),
            ),
            (false, false, true, false) => IntervalNumber::new(
                (self.low * rhs.high).next_up(),
                (self.high * rhs.high).next_up(),
            ),
            (false, false, false, true) => IntervalNumber::new(
                (self.high * rhs.low).next_up(),
                (self.low * rhs.low).next_up(),
            ),
            (false, false, false, false) => IntervalNumber::new(
                ((self.low * rhs.high).max(self.high * rhs.low)).next_up(),
                ((self.low * rhs.low).max(self.high * rhs.high)).next_up(),
            ),
            _ => IntervalNumber::default(),
        }
    }
}
