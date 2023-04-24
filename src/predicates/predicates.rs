#![allow(non_snake_case)]

use bumpalo::{collections::Vec, Bump};

struct Bound {
    epsilon: f64,
    splitter: f64,
    resulterrbound: f64,
    ccwerrboundA: f64,
    ccwerrboundB: f64,
    ccwerrboundC: f64,
    o3derrboundA: f64,
    o3derrboundB: f64,
    o3derrboundC: f64,
    iccerrboundA: f64,
    iccerrboundB: f64,
    iccerrboundC: f64,
    isperrboundA: f64,
    isperrboundB: f64,
    isperrboundC: f64,
}

const fn exact_init() -> Bound {
    let mut every_other = true;
    let half = 0.5;
    let mut epsilon = 1.0;
    let mut splitter = 1.0;
    let mut check = 1.0;
    /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    /*   one without causing round off.  (Also check if the sum is equal to  */
    /*   the previous sum, for machines that round up instead of using exact */
    /*   rounding.  Not that this library will work on such machines anyway. */
    loop {
        let lastcheck = check;
        epsilon *= half;

        if every_other {
            splitter *= 2.0;
        }

        every_other = !every_other;
        check = 1.0 + epsilon;
        if check != 1.0 && check != lastcheck {
            break;
        }
    }
    splitter += 1.0;

    /* Error bounds for orientation and incircle tests. */
    let resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
    let ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
    let ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
    let ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    let o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
    let o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
    let o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
    let iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
    let iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
    let iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
    let isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
    let isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
    let isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
    return Bound {
        epsilon,
        splitter,
        resulterrbound,
        ccwerrboundA,
        ccwerrboundB,
        ccwerrboundC,
        o3derrboundA,
        o3derrboundB,
        o3derrboundC,
        iccerrboundA,
        iccerrboundB,
        iccerrboundC,
        isperrboundA,
        isperrboundB,
        isperrboundC,
    };
}

const B: Bound = exact_init();

#[inline(always)]
pub const fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bvirt = x - a;
    let y = b - bvirt;
    (x, y)
}

#[inline(always)]
pub const fn fast_two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    let bvirt = a - x;
    let y = bvirt - b;
    (x, y)
}

#[inline(always)]
pub const fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bvirt = x - a;
    let avirt = x - bvirt;
    let bround = b - bvirt;
    let around = a - avirt;
    let y = around + bround;
    (x, y)
}

#[inline(always)]
pub const fn two_diff_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt = a - x;
    let avirt = x + bvirt;
    let bround = bvirt - b;
    let around = a - avirt;
    around + bround
}

#[inline(always)]
pub const fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    let bvirt = a - x;
    let avirt = x + bvirt;
    let bround = bvirt - b;
    let around = a - avirt;
    (x, around + bround)
}

#[inline(always)]
const fn split(a: f64) -> (f64, f64) {
    let c = B.splitter * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;
    (ahi, alo)
}

#[inline(always)]
pub const fn two_product(a: f64, b: f64) -> (f64, f64) {
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let x = a * b;
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    let y = (alo * blo) - err3;
    (x, y)
}

#[inline(always)]
pub const fn two_product_pre_split(a: f64, b: f64, bhi: f64, blo: f64) -> (f64, f64) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);
    let y = (alo * blo) - err3;
    (x, y)
}

#[inline(always)]
const fn two_one_sum(a1: f64, a0: f64, b: f64) -> (f64, f64, f64) {
    let (_i, x0) = two_sum(a0, b);
    let (x2, x1) = two_sum(a1, _i);
    (x2, x1, x0)
}

#[inline(always)]
pub const fn two_two_sum(a1: f64, a0: f64, b1: f64, b0: f64) -> (f64, f64, f64, f64) {
    let (_j, _0, x0) = two_one_sum(a1, a0, b0);
    let (x3, x2, x1) = two_one_sum(_j, _0, b1);
    (x3, x2, x1, x0)
}

#[inline(always)]
const fn two_one_diff(a1: f64, a0: f64, b: f64) -> (f64, f64, f64) {
    let (_i, x0) = two_diff(a0, b);
    let (x2, x1) = two_sum(a1, _i);
    (x2, x1, x0)
}

#[inline(always)]
pub const fn two_two_diff(a1: f64, a0: f64, b1: f64, b0: f64) -> (f64, f64, f64, f64) {
    let (_j, _0, x0) = two_one_diff(a1, a0, b0);
    let (x3, x2, x1) = two_one_diff(_j, _0, b1);
    (x3, x2, x1, x0)
}

#[inline]
pub fn grow_expansion<'a>(earr: &[f64], b: f64, bump: &'a Bump) -> Vec<'a, f64> {
    let mut q = b;
    let mut harr = Vec::from_iter_in(
        earr.iter().map(|&enow| {
            let (qnew, h) = two_sum(q, enow);
            q = qnew;
            h
        }),
        bump,
    );
    harr.push(q);
    harr
}

#[inline]
pub fn grow_expansion_zeroelim<'a>(earr: &[f64], b: f64, bump: &'a Bump) -> Vec<'a, f64> {
    let mut q = b;
    let mut h = Vec::from_iter_in(
        earr.iter().filter_map(|&enow| {
            let (qnew, r) = two_sum(q, enow);
            q = qnew;
            if r != 0.0 {
                Some(r)
            } else {
                None
            }
        }),
        bump,
    );
    if q != 0.0 {
        h.push(q);
    }
    h
}

#[inline]
pub fn expansion_sum<'a>(earr: &[f64], farr: &[f64], bump: &'a Bump) -> Vec<'a, f64> {
    let mut q = farr[0];
    let mut h_arr = Vec::with_capacity_in(earr.len() + farr.len(), bump);
    h_arr.extend(earr.iter().map(|&enow| {
        let (qnew, r) = two_sum(q, enow);
        q = qnew;
        r
    }));
    h_arr.push(q);
    earr.iter().enumerate().skip(1).for_each(|(start, &e)| {
        q = e;
        for h in &mut h_arr[start..] {
            let (qnew, enew) = two_sum(q, *h);
            q = qnew;
            *h = enew;
        }
        h_arr.push(q);
    });
    h_arr
}

struct TwoArr<'a> {
    arr1: &'a [f64],
    arr2: &'a [f64],
}

impl<'a> TwoArr<'a> {
    #[inline(always)]
    fn first(&self) -> f64 {
        self.arr1[0]
    }

    #[inline(always)]
    fn second(&self) -> f64 {
        self.arr2[0]
    }

    #[inline(always)]
    fn first_advance(&mut self) {
        self.arr1 = &self.arr1[1..];
    }

    #[inline(always)]
    fn second_advance(&mut self) {
        self.arr2 = &self.arr2[1..];
    }
}

impl<'a> Iterator for TwoArr<'a> {
    type Item = f64;
    fn next(&mut self) -> Option<Self::Item> {
        let item: f64;
        if !self.arr1.is_empty() && !self.arr2.is_empty() {
            let (first, second) = (self.first(), self.second());
            if (first > second) == (first > -second) {
                item = second;
                self.second_advance();
            } else {
                item = first;
                self.first_advance();
            }
        } else {
            if !self.arr1.is_empty() {
                item = self.first();
                self.first_advance();
            } else if !self.arr2.is_empty() {
                item = self.second();
            } else {
                return None;
            }
        }
        Some(item)
    }
}

pub fn fast_expansion_sum_zeroelim<'a>(earr: &[f64], farr: &[f64], bump: &'a Bump) -> Vec<'a, f64> {
    let mut iter = TwoArr {
        arr1: earr,
        arr2: farr,
    };

    let mut q = iter.next().unwrap();
    if let Some(now) = iter.next() {
        // init
        let (qnew, h) = fast_two_sum(now, q);
        q = qnew;
        let mut harr = Vec::with_capacity_in(earr.len() * farr.len(), bump);
        if h != 0.0 {
            harr.push(h);
        }

        for now in iter {
            let (qnew, h) = two_sum(q, now);
            q = qnew;
            if h != 0.0 {
                harr.push(h);
            }
        }
        if q != 0.0 || harr.is_empty() {
            harr.push(q)
        }
        harr
    } else {
        Vec::new_in(bump)
    }
}

#[inline]
pub fn fast_expansion_diff_zeroelim<'a>(
    earr: &[f64],
    farr: &[f64],
    bump: &'a Bump,
) -> Vec<'a, f64> {
    let farr_oppo = Vec::from_iter_in(farr.iter().map(|&f| -f), bump);
    fast_expansion_sum_zeroelim(earr, &farr_oppo, bump)
}

pub fn scale_expansion_zeroelim<'a>(earr: &[f64], b: f64, bump: &'a Bump) -> Vec<'a, f64> {
    let (bhi, blo) = split(b);
    let (&efirst, earr) = earr.split_first().unwrap();
    let (mut q, h) = two_product_pre_split(efirst, b, bhi, blo);
    let mut harr = Vec::new_in(bump);
    if h != 0.0 {
        harr.push(h);
    }

    harr.extend(
        earr.iter()
            .map(|&e| {
                let (p1, p0) = two_product_pre_split(e, b, bhi, blo);
                let (sum, mut h) = two_sum(q, p0);
                let mut ret = Vec::with_capacity_in(2, bump);
                if h != 0.0 {
                    ret.push(h);
                }

                (q, h) = fast_two_sum(p1, sum);
                if h != 0.0 {
                    ret.push(h);
                }
                ret
            })
            .flatten(),
    );
    if q != 0.0 || harr.is_empty() {
        harr.push(q);
    }
    harr
}

#[inline(always)]
pub fn estimate(arr: &[f64]) -> f64 {
    arr.iter().sum()
}

fn orient2d_adapt<'a>(pa: &[f64], pb: &[f64], pc: &[f64], detsum: f64, bump: &'a Bump) -> f64 {
    let acx = pa[0] - pc[0];
    let bcx = pb[0] - pc[0];
    let acy = pa[1] - pc[1];
    let bcy = pb[1] - pc[1];

    let (detleft, detlefttail) = two_product(acx, bcy);
    let (detright, detrighttail) = two_product(acy, bcx);
    let mut b = bumpalo::vec![in bump; 0.0; 4];
    (b[3], b[2], b[1], b[0]) = two_two_diff(detleft, detlefttail, detright, detrighttail);

    let mut det = estimate(&b);
    let errbound = B.ccwerrboundB * detsum;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let acxtail = two_diff_tail(pa[0], pc[0], acx);
    let bcxtail = two_diff_tail(pb[0], pc[0], bcx);
    let acytail = two_diff_tail(pa[1], pc[1], acy);
    let bcytail = two_diff_tail(pb[1], pc[1], bcy);

    if (acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) && (bcytail == 0.0) {
        return det;
    }

    let errbound = B.ccwerrboundC * detsum + B.resulterrbound * det.abs();
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let (s1, s0) = two_product(acxtail, bcy);
    let (t1, t0) = two_product(acytail, bcx);
    let mut u = bumpalo::vec![in bump; 0.0; 4];
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let c1 = fast_expansion_sum_zeroelim(&b, &u, bump);

    let (s1, s0) = two_product(acx, bcytail);
    let (t1, t0) = two_product(acy, bcxtail);
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let c2 = fast_expansion_sum_zeroelim(&c1, &u, bump);

    let (s1, s0) = two_product(acxtail, bcytail);
    let (t1, t0) = two_product(acytail, bcxtail);
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let d = fast_expansion_sum_zeroelim(&c2, &u, bump);

    *d.last().unwrap()
}

pub fn orien2d<'a>(pa: &[f64], pb: &[f64], pc: &[f64], bump: &'a Bump) -> f64 {
    let detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    let detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    let det = detleft - detright;
    let detsum: f64;

    if detleft > 0.0 {
        if detright <= 0.0 {
            return det;
        } else {
            detsum = detleft + detright;
        }
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return det;
        } else {
            detsum = -detleft - detright;
        }
    } else {
        return det;
    }

    let errbound = B.ccwerrboundA * detsum;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    orient2d_adapt(pa, pb, pc, detsum, bump)
}
