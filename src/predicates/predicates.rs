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
pub const fn two_one_product(a1: f64, a0: f64, b: f64) -> (f64, f64, f64, f64) {
    let (bhi, blo) = split(b);
    let (_i, x0) = two_product_pre_split(a0, b, bhi, blo);
    let (_j, _0) = two_product_pre_split(a1, b, bhi, blo);
    let (_k, x1) = two_sum(_i, _0);
    let (x3, x2) = fast_two_sum(_j, _k);
    (x3, x2, x1, x0)
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

fn orient3d_adapt<'a>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    permanent: f64,
    bump: &'a Bump,
) -> f64 {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdz = pb[2] - pd[2];
    let cdz = pc[2] - pd[2];

    let (bdxcdy1, bdxcdy0) = two_product(bdx, cdy);
    let (cdxbdy1, cdxbdy0) = two_product(cdx, bdy);
    let mut bc = bumpalo::vec![in bump; 0.0; 4];
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let adet = scale_expansion_zeroelim(&bc, adz, bump);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let mut ca = bumpalo::vec![in bump; 0.0; 4];
    (ca[3], ca[2], ca[1], ca[0]) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let bdet = scale_expansion_zeroelim(&ca, bdz, bump);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let mut ab = bumpalo::vec![in bump; 0.0; 4];
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let cdet = scale_expansion_zeroelim(&ab, cdz, bump);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, bump);
    let mut fin = fast_expansion_sum_zeroelim(&abdet, &cdet, bump);

    let mut det = estimate(&fin);
    let errbound = B.o3derrboundB * permanent;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let adxtail = two_diff_tail(pa[0], pd[0], adx);
    let bdxtail = two_diff_tail(pb[0], pd[0], bdx);
    let cdxtail = two_diff_tail(pc[0], pd[0], cdx);
    let adytail = two_diff_tail(pa[1], pd[1], ady);
    let bdytail = two_diff_tail(pb[1], pd[1], bdy);
    let cdytail = two_diff_tail(pc[1], pd[1], cdy);
    let adztail = two_diff_tail(pa[2], pd[2], adz);
    let bdztail = two_diff_tail(pb[2], pd[2], bdz);
    let cdztail = two_diff_tail(pc[2], pd[2], cdz);

    if (adxtail == 0.0)
        && (bdxtail == 0.0)
        && (cdxtail == 0.0)
        && (adytail == 0.0)
        && (bdytail == 0.0)
        && (cdytail == 0.0)
        && (adztail == 0.0)
        && (bdztail == 0.0)
        && (cdztail == 0.0)
    {
        return det;
    }

    let errbound = B.o3derrboundC * permanent + B.resulterrbound * det.abs();
    det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + adztail * (bdx * cdy - bdy * cdx))
        + (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + bdztail * (cdx * ady - cdy * adx))
        + (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + cdztail * (adx * bdy - ady * bdx));
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    fn helper(
        adxtail: f64,
        adytail: f64,
        bdx: f64,
        bdy: f64,
        cdx: f64,
        cdy: f64,
        bump: &Bump,
    ) -> (Vec<f64>, Vec<f64>) {
        if adxtail == 0.0 {
            if adytail == 0.0 {
                (bumpalo::vec![in bump; 0.0], bumpalo::vec![in bump; 0.0])
            } else {
                let (mut at_b, mut at_c) = (
                    bumpalo::vec![in bump; 0.0, 0.0],
                    bumpalo::vec![in bump; 0.0, 0.0],
                );
                (at_b[1], at_b[0]) = two_product(-adytail, bdx);
                (at_c[1], at_c[0]) = two_product(adytail, cdx);
                (at_b, at_c)
            }
        } else {
            if adytail == 0.0 {
                let (mut at_b, mut at_c) = (
                    bumpalo::vec![in bump; 0.0, 0.0],
                    bumpalo::vec![in bump; 0.0, 0.0],
                );
                (at_b[1], at_b[0]) = two_product(adxtail, bdy);
                (at_c[1], at_c[0]) = two_product(-adxtail, cdy);
                (at_b, at_c)
            } else {
                let (adxt_bdy1, adxt_bdy0) = two_product(adxtail, bdy);
                let (adyt_bdx1, adyt_bdx0) = two_product(adytail, bdx);
                let (mut at_b, mut at_c) = (
                    bumpalo::vec![in bump; 0.0; 4],
                    bumpalo::vec![in bump; 0.0; 4],
                );
                (at_b[3], at_b[2], at_b[1], at_b[0]) =
                    two_two_diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0);
                let (adyt_cdx1, adyt_cdx0) = two_product(adytail, cdx);
                let (adxt_cdy1, adxt_cdy0) = two_product(adxtail, cdy);
                (at_c[3], at_c[2], at_c[1], at_c[0]) =
                    two_two_diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0);
                (at_b, at_c)
            }
        }
    }

    let (at_b, at_c) = helper(adxtail, adytail, bdx, bdy, cdx, cdy, bump);
    let (bt_c, bt_a) = helper(bdxtail, bdytail, cdx, cdy, adx, ady, bump);
    let (ct_a, ct_b) = helper(cdxtail, cdytail, adx, ady, bdx, bdy, bump);

    let xyt = Vec::from_iter_in(
        [(bt_c, ct_b, adz), (ct_a, at_c, bdz), (at_b, bt_a, cdz)]
            .into_iter()
            .map(|tuple| {
                let sum = fast_expansion_sum_zeroelim(&tuple.0, &tuple.1, bump);
                let w = scale_expansion_zeroelim(&sum, tuple.2, bump);
                fin = fast_expansion_sum_zeroelim(&fin, &w, bump);
                sum
            }),
        bump,
    );

    if adztail != 0.0 {
        let v = scale_expansion_zeroelim(&bc, adztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }

    if bdztail != 0.0 {
        let v = scale_expansion_zeroelim(&ca, bdztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }

    if cdztail != 0.0 {
        let v = scale_expansion_zeroelim(&ab, cdztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }

    fn helper2<'a>(
        adxtail: f64,
        bdytail: f64,
        bdztail: f64,
        bdz: f64,
        cdz: f64,
        cdytail: f64,
        cdztail: f64,
        fin: &mut Vec<'a, f64>,
        bump: &'a Bump,
    ) {
        if adxtail != 0.0 {
            if bdytail != 0.0 {
                let (adxt_bdyt1, adxt_bdyt0) = two_product(adxtail, bdytail);
                let mut u = bumpalo::vec![in bump; 0.0; 4];
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdz);
                *fin = fast_expansion_sum_zeroelim(fin, &u, bump);
                if cdztail != 0.0 {
                    (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail);
                    *fin = fast_expansion_sum_zeroelim(fin, &u, bump);
                }
            }
            if cdytail != 0.0 {
                let (adxt_cdyt1, adxt_cdyt0) = two_product(-adxtail, cdytail);
                let mut u = bumpalo::vec![in bump; 0.0; 4];
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdz);
                *fin = fast_expansion_sum_zeroelim(fin, &u, bump);
                if bdztail != 0.0 {
                    (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail);
                    *fin = fast_expansion_sum_zeroelim(fin, &u, bump);
                }
            }
        }
    }

    helper2(
        adxtail, bdytail, bdztail, bdz, cdz, cdytail, cdztail, &mut fin, bump,
    );
    helper2(
        bdxtail, cdytail, cdztail, cdz, adz, adytail, adztail, &mut fin, bump,
    );
    helper2(
        cdztail, adytail, adztail, adz, bdz, bdytail, bdztail, &mut fin, bump,
    );

    if adztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[0], adztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }

    if bdztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[1], bdztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }

    if cdztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[2], cdztail, bump);
        fin = fast_expansion_sum_zeroelim(&fin, &v, bump);
    }
    *fin.last().unwrap()
}
