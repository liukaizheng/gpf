use std::ops::{
    Add, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, Index, Not, Rem, Shl, Shr, Sub,
};

use tinyvec::{Array, TinyVec};

pub trait BitBlock:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + Rem<Output = Self>
    + BitOrAssign
    + BitAndAssign
    + PartialEq
    + Eq
    + std::hash::Hash
{
    /// How many bits it has
    fn bits() -> usize;
    /// How many bytes it has
    #[inline]
    fn bytes() -> usize {
        Self::bits() / 8
    }
    /// Convert a byte into this type (lowest-order bits set)
    fn from_byte(byte: u8) -> Self;
    /// Count the number of 1's in the bitwise repr
    fn count_ones(self) -> usize;
    /// Count the number of 0's in the bitwise repr
    fn count_zeros(self) -> usize {
        Self::bits() - self.count_ones()
    }
    /// Get `0`
    fn zero() -> Self;
    /// Get `1`
    fn one() -> Self;
}

macro_rules! bit_block_impl {
    ($(($t: ident, $size: expr)),*) => ($(
        impl BitBlock for $t {
            #[inline]
            fn bits() -> usize { $size }
            #[inline]
            fn from_byte(byte: u8) -> Self { $t::from(byte) }
            #[inline]
            fn count_ones(self) -> usize { self.count_ones() as usize }
            #[inline]
            fn count_zeros(self) -> usize { self.count_zeros() as usize }
            #[inline]
            fn one() -> Self { 1 }
            #[inline]
            fn zero() -> Self { 0 }
        }
    )*)
}

bit_block_impl! {
    (u8, 8),
    (u16, 16),
    (u32, 32),
    (u64, 64),
    (usize, core::mem::size_of::<usize>() * 8)
}

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct Bitmask<A: Array = [usize; 1]>
where
    A::Item: BitBlock,
{
    data: TinyVec<A>,
}

impl<A: Array> Bitmask<A>
where
    A::Item: BitBlock,
{
    pub fn new(len: usize) -> Self {
        let n = (len + A::Item::bits() - 1) / A::Item::bits();

        let data = TinyVec::from_iter((0..n).map(|_| A::Item::zero()));
        Self { data }
    }

    #[inline]
    pub fn set(&mut self, i: usize) {
        let block = i / A::Item::bits();
        let bit = i % A::Item::bits();
        self.data[block] |= A::Item::one() << bit;
    }

    #[inline]
    pub fn clear(&mut self, i: usize) {
        let block = i / A::Item::bits();
        let bit = i % A::Item::bits();
        self.data[block] &= !(A::Item::one() << bit);
    }

    #[inline]
    pub fn n_elements(&self) -> usize {
        self.data.iter().map(|&x| x.count_ones()).sum()
    }

    #[inline]
    pub fn n_common(&self, other: &Self) -> usize {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| (a & b).count_ones())
            .sum()
    }
}

impl<A: Array> Index<usize> for Bitmask<A>
where
    A::Item: BitBlock,
{
    type Output = bool;

    fn index(&self, index: usize) -> &Self::Output {
        let block = index / A::Item::bits();
        let bit = index % A::Item::bits();
        let mask = A::Item::one() << bit;
        let data = self.data[block];
        if (data & mask) == A::Item::zero() {
            &false
        } else {
            &true
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitmask() {
        let mut bm1 = Bitmask::<[u8; 1]>::new(10);
        assert_eq!(bm1.data.len(), 2);
        [1, 2, 7, 9, 7].iter().for_each(|&i| bm1.set(i));
        assert_eq!(bm1[0], false);
        assert_eq!(bm1[2], true);

        let mut bm2 = Bitmask::<[u8; 1]>::new(10);
        [2, 7, 9].iter().for_each(|&i| bm2.set(i));

        assert_eq!(bm1.n_common(&bm2), 3);
    }
}
