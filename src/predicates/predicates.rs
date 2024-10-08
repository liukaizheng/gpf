use std::alloc::Allocator;

struct Bound {
    splitter: f64,
    result_err_bound: f64,
    ccw_err_bounda: f64,
    ccw_err_boundb: f64,
    ccw_err_boundc: f64,
    o3d_err_bounda: f64,
    o3d_err_boundb: f64,
    o3d_err_boundc: f64,
    icc_err_bounda: f64,
    icc_err_boundb: f64,
    icc_err_boundc: f64,
    isp_err_bounda: f64,
    isp_err_boundb: f64,
    isp_err_boundc: f64,
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
        let last_check = check;
        epsilon *= half;

        if every_other {
            splitter *= 2.0;
        }

        every_other = !every_other;
        check = 1.0 + epsilon;
        if check == 1.0 || check == last_check {
            break;
        }
    }
    splitter += 1.0;

    /* Error bounds for orientation and incircle tests. */
    let result_err_bound = (3.0 + 8.0 * epsilon) * epsilon;
    let ccw_err_bounda = (3.0 + 16.0 * epsilon) * epsilon;
    let ccw_err_boundb = (2.0 + 12.0 * epsilon) * epsilon;
    let ccw_err_boundc = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    let o3d_err_bounda = (7.0 + 56.0 * epsilon) * epsilon;
    let o3d_err_boundb = (3.0 + 28.0 * epsilon) * epsilon;
    let o3d_err_boundc = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
    let icc_err_bounda = (10.0 + 96.0 * epsilon) * epsilon;
    let icc_err_boundb = (4.0 + 48.0 * epsilon) * epsilon;
    let icc_err_boundc = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
    let isp_err_bounda = (16.0 + 224.0 * epsilon) * epsilon;
    let isp_err_boundb = (5.0 + 72.0 * epsilon) * epsilon;
    let isp_err_boundc = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
    Bound {
        splitter,
        result_err_bound,
        ccw_err_bounda,
        ccw_err_boundb,
        ccw_err_boundc,
        o3d_err_bounda,
        o3d_err_boundb,
        o3d_err_boundc,
        icc_err_bounda,
        icc_err_boundb,
        icc_err_boundc,
        isp_err_bounda,
        isp_err_boundb,
        isp_err_boundc,
    }
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
pub const fn square(a: f64) -> (f64, f64) {
    let x = a * a;
    let (ahi, alo) = split(a);
    let err1 = x - (ahi * ahi);
    let err3 = err1 - ((ahi + ahi) * alo);
    let y = (alo * alo) - err3;
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
pub fn grow_expansion<A: Allocator>(earr: &[f64], b: f64, allocator: A) -> Vec<f64, A> {
    let mut q = b;
    let mut harr = Vec::new_in(allocator);
    harr.extend(earr.iter().map(|&enow| {
        let (qnew, h) = two_sum(q, enow);
        q = qnew;
        h
    }));
    harr.push(q);
    harr
}

#[inline]
pub fn grow_expansion_zeroelim<A: Allocator>(earr: &[f64], b: f64, allocator: A) -> Vec<f64, A> {
    let mut q = b;
    let mut h = Vec::new_in(allocator);
    h.extend(earr.iter().filter_map(|&enow| {
        let (qnew, r) = two_sum(q, enow);
        q = qnew;
        if r != 0.0 {
            Some(r)
        } else {
            None
        }
    }));
    if q != 0.0 {
        h.push(q);
    }
    h
}

#[inline]
pub fn expansion_sum<A: Allocator>(earr: &[f64], farr: &[f64], allocator: A) -> Vec<f64, A> {
    let mut q = farr[0];
    let mut h_arr = Vec::with_capacity_in(earr.len() + farr.len(), allocator);
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
                self.second_advance();
            } else {
                return None;
            }
        }
        Some(item)
    }
}

pub fn fast_expansion_sum_zeroelim<A: Allocator>(
    earr: &[f64],
    farr: &[f64],
    allocator: A,
) -> Vec<f64, A> {
    let mut iter = TwoArr {
        arr1: earr,
        arr2: farr,
    };

    let mut q = iter.next().unwrap();
    if let Some(now) = iter.next() {
        // init
        let (qnew, h) = fast_two_sum(now, q);
        q = qnew;
        let mut harr = Vec::with_capacity_in(earr.len() + farr.len(), allocator);
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
        Vec::new_in(allocator)
    }
}

#[inline]
pub fn fast_expansion_diff_zeroelim<A: Allocator + Copy>(
    earr: &[f64],
    farr: &[f64],
    allocator: A,
) -> Vec<f64, A> {
    let mut farr_oppo = Vec::new_in(allocator);
    farr_oppo.extend(farr.iter().map(|&f| -f));
    fast_expansion_sum_zeroelim(earr, &farr_oppo, allocator)
}

pub fn scale_expansion_zeroelim<A: Allocator + Copy>(
    earr: &[f64],
    b: f64,
    allocator: A,
) -> Vec<f64, A> {
    let (bhi, blo) = split(b);
    let (&efirst, earr) = earr.split_first().unwrap();
    let (mut q, h) = two_product_pre_split(efirst, b, bhi, blo);
    let mut harr = Vec::new_in(allocator);
    if h != 0.0 {
        harr.push(h);
    }

    harr.extend(
        earr.iter()
            .map(|&e| {
                let (p1, p0) = two_product_pre_split(e, b, bhi, blo);
                let (sum, mut h) = two_sum(q, p0);
                let mut ret = Vec::with_capacity_in(2, allocator);
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
fn mul_expansion_zeroelim_imp<A: Allocator + Copy>(
    earr: &[f64],
    farr: &[f64],
    allocator: A,
) -> Vec<f64, A> {
    earr.iter()
        .map(|&s| scale_expansion_zeroelim(&farr, s, allocator))
        .reduce(|a, b| fast_expansion_sum_zeroelim(&a, &b, allocator))
        .unwrap()
}

#[inline(always)]
pub fn mul_expansion_zeroelim<A: Allocator + Copy>(
    earr: &[f64],
    farr: &[f64],
    allocator: A,
) -> Vec<f64, A> {
    if earr.len() == 1 {
        scale_expansion_zeroelim(farr, earr[0], allocator)
    } else if farr.len() == 1 {
        scale_expansion_zeroelim(earr, farr[0], allocator)
    } else {
        if earr.len() < farr.len() {
            mul_expansion_zeroelim_imp(earr, farr, allocator)
        } else {
            mul_expansion_zeroelim_imp(farr, earr, allocator)
        }
    }
}

#[inline(always)]
pub fn expansion_invert(arr: &mut [f64]) {
    for a in arr {
        *a = -*a;
    }
}

#[inline(always)]
pub fn estimate(arr: &[f64]) -> f64 {
    arr.iter().sum()
}

fn orient2d_adapt<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    detsum: f64,
    allocator: A,
) -> f64 {
    let acx = pa[0] - pc[0];
    let bcx = pb[0] - pc[0];
    let acy = pa[1] - pc[1];
    let bcy = pb[1] - pc[1];

    let (detleft, detlefttail) = two_product(acx, bcy);
    let (detright, detrighttail) = two_product(acy, bcx);
    let mut b = std::vec::from_elem_in(0.0, 4, allocator);
    (b[3], b[2], b[1], b[0]) = two_two_diff(detleft, detlefttail, detright, detrighttail);

    let mut det = estimate(&b);
    let errbound = B.ccw_err_boundb * detsum;
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

    let errbound = B.ccw_err_boundc * detsum + B.result_err_bound * det.abs();
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let (s1, s0) = two_product(acxtail, bcy);
    let (t1, t0) = two_product(acytail, bcx);
    let mut u = std::vec::from_elem_in(0.0, 4, allocator);
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let c1 = fast_expansion_sum_zeroelim(&b, &u, allocator);

    let (s1, s0) = two_product(acx, bcytail);
    let (t1, t0) = two_product(acy, bcxtail);
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let c2 = fast_expansion_sum_zeroelim(&c1, &u, allocator);

    let (s1, s0) = two_product(acxtail, bcytail);
    let (t1, t0) = two_product(acytail, bcxtail);
    (u[3], u[2], u[1], u[0]) = two_two_diff(s1, s0, t1, t0);
    let d = fast_expansion_sum_zeroelim(&c2, &u, allocator);

    *d.last().unwrap()
}

pub fn orient2d<A: Allocator + Copy>(pa: &[f64], pb: &[f64], pc: &[f64], allocator: A) -> f64 {
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

    let errbound = B.ccw_err_bounda * detsum;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    orient2d_adapt(pa, pb, pc, detsum, allocator)
}

fn orient3d_adapt<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    permanent: f64,
    allocator: A,
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
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let adet = scale_expansion_zeroelim(&bc, adz, allocator);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let mut ca = std::vec::from_elem_in(0.0, 4, allocator);
    (ca[3], ca[2], ca[1], ca[0]) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let bdet = scale_expansion_zeroelim(&ca, bdz, allocator);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let cdet = scale_expansion_zeroelim(&ab, cdz, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let mut fin = fast_expansion_sum_zeroelim(&abdet, &cdet, allocator);

    let mut det = estimate(&fin);
    let errbound = B.o3d_err_boundb * permanent;
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

    let errbound = B.o3d_err_boundc * permanent + B.result_err_bound * det.abs();
    det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + adztail * (bdx * cdy - bdy * cdx))
        + (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + bdztail * (cdx * ady - cdy * adx))
        + (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + cdztail * (adx * bdy - ady * bdx));
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    fn helper1<A: Allocator + Copy>(
        adxtail: f64,
        adytail: f64,
        bdx: f64,
        bdy: f64,
        cdx: f64,
        cdy: f64,
        allocator: A,
    ) -> (Vec<f64, A>, Vec<f64, A>) {
        if adxtail == 0.0 {
            if adytail == 0.0 {
                (
                    std::vec::from_elem_in(0.0, 4, allocator),
                    std::vec::from_elem_in(0.0, 4, allocator),
                )
            } else {
                let (mut at_b, mut at_c) = (
                    std::vec::from_elem_in(0.0, 4, allocator),
                    std::vec::from_elem_in(0.0, 4, allocator),
                );
                (at_b[1], at_b[0]) = two_product(-adytail, bdx);
                (at_c[1], at_c[0]) = two_product(adytail, cdx);
                (at_b, at_c)
            }
        } else {
            if adytail == 0.0 {
                let (mut at_b, mut at_c) = (
                    std::vec::from_elem_in(0.0, 4, allocator),
                    std::vec::from_elem_in(0.0, 4, allocator),
                );
                (at_b[1], at_b[0]) = two_product(adxtail, bdy);
                (at_c[1], at_c[0]) = two_product(-adxtail, cdy);
                (at_b, at_c)
            } else {
                let (adxt_bdy1, adxt_bdy0) = two_product(adxtail, bdy);
                let (adyt_bdx1, adyt_bdx0) = two_product(adytail, bdx);
                let (mut at_b, mut at_c) = (
                    std::vec::from_elem_in(0.0, 4, allocator),
                    std::vec::from_elem_in(0.0, 4, allocator),
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

    let (at_b, at_c) = helper1(adxtail, adytail, bdx, bdy, cdx, cdy, allocator);
    let (bt_c, bt_a) = helper1(bdxtail, bdytail, cdx, cdy, adx, ady, allocator);
    let (ct_a, ct_b) = helper1(cdxtail, cdytail, adx, ady, bdx, bdy, allocator);

    let mut xyt = Vec::new_in(allocator);
    xyt.extend(
        [(bt_c, ct_b, adz), (ct_a, at_c, bdz), (at_b, bt_a, cdz)]
            .into_iter()
            .map(|tuple| {
                let sum = fast_expansion_sum_zeroelim(&tuple.0, &tuple.1, allocator);
                let w = scale_expansion_zeroelim(&sum, tuple.2, allocator);
                fin = fast_expansion_sum_zeroelim(&fin, &w, allocator);
                sum
            }),
    );

    if adztail != 0.0 {
        let v = scale_expansion_zeroelim(&bc, adztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }

    if bdztail != 0.0 {
        let v = scale_expansion_zeroelim(&ca, bdztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }

    if cdztail != 0.0 {
        let v = scale_expansion_zeroelim(&ab, cdztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }

    fn helper2<A: Allocator + Copy>(
        adxtail: f64,
        bdytail: f64,
        bdztail: f64,
        bdz: f64,
        cdz: f64,
        cdytail: f64,
        cdztail: f64,
        fin: &mut Vec<f64, A>,
        allocator: A,
    ) {
        if adxtail != 0.0 {
            if bdytail != 0.0 {
                let (adxt_bdyt1, adxt_bdyt0) = two_product(adxtail, bdytail);
                let mut u = std::vec::from_elem_in(0.0, 4, allocator);
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdz);
                *fin = fast_expansion_sum_zeroelim(fin, &u, allocator);
                if cdztail != 0.0 {
                    (u[3], u[2], u[1], u[0]) = two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail);
                    *fin = fast_expansion_sum_zeroelim(fin, &u, allocator);
                }
            }
            if cdytail != 0.0 {
                let (adxt_cdyt1, adxt_cdyt0) = two_product(-adxtail, cdytail);
                let mut u = std::vec::from_elem_in(0.0, 4, allocator);
                (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdz);
                *fin = fast_expansion_sum_zeroelim(fin, &u, allocator);
                if bdztail != 0.0 {
                    (u[3], u[2], u[1], u[0]) = two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail);
                    *fin = fast_expansion_sum_zeroelim(fin, &u, allocator);
                }
            }
        }
    }

    helper2(
        adxtail, bdytail, bdztail, bdz, cdz, cdytail, cdztail, &mut fin, allocator,
    );
    helper2(
        bdxtail, cdytail, cdztail, cdz, adz, adytail, adztail, &mut fin, allocator,
    );
    helper2(
        cdztail, adytail, adztail, adz, bdz, bdytail, bdztail, &mut fin, allocator,
    );

    if adztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[0], adztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }

    if bdztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[1], bdztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }

    if cdztail != 0.0 {
        let v = scale_expansion_zeroelim(&xyt[2], cdztail, allocator);
        fin = fast_expansion_sum_zeroelim(&fin, &v, allocator);
    }
    *fin.last().unwrap()
}

pub fn orient3d<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    allocator: A,
) -> f64 {
    let adx = pa[0] - pd[0];
    let ady = pa[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdx = pb[0] - pd[0];
    let bdy = pb[1] - pd[1];
    let bdz = pb[2] - pd[2];
    let cdx = pc[0] - pd[0];
    let cdy = pc[1] - pd[1];
    let cdz = pc[2] - pd[2];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);

    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
        + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
        + (adxbdy.abs() + bdxady.abs()) * cdz.abs();
    let err_bound = B.o3d_err_bounda * permanent;
    if (det > err_bound) || (-det > err_bound) {
        return det;
    }

    orient3d_adapt(pa, pb, pc, pd, permanent, allocator)
}

pub fn orient3d_fast(pa: &[f64], pb: &[f64], pc: &[f64], pd: &[f64]) -> f64 {
    let adx = pa[0] - pd[0];
    let ady = pa[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdx = pb[0] - pd[0];
    let bdy = pb[1] - pd[1];
    let bdz = pb[2] - pd[2];
    let cdx = pc[0] - pd[0];
    let cdy = pc[1] - pd[1];
    let cdz = pc[2] - pd[2];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady)
}

fn incircle_adapt<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    permanent: f64,
    allocator: A,
) -> f64 {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];

    let (bdxcdy1, bdxcdy0) = two_product(bdx, cdy);
    let (cdxbdy1, cdxbdy0) = two_product(cdx, bdy);
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let axbc = scale_expansion_zeroelim(&bc, adx, allocator);
    let axxbc = scale_expansion_zeroelim(&axbc, adx, allocator);
    let aybc = scale_expansion_zeroelim(&bc, ady, allocator);
    let ayybc = scale_expansion_zeroelim(&aybc, ady, allocator);
    let adet = fast_expansion_sum_zeroelim(&axxbc, &ayybc, allocator);

    let (cdxady1, cdxady0) = two_product(cdx, ady);
    let (adxcdy1, adxcdy0) = two_product(adx, cdy);
    let mut ca = std::vec::from_elem_in(0.0, 4, allocator);
    (ca[3], ca[2], ca[1], ca[0]) = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let bxca = scale_expansion_zeroelim(&ca, bdx, allocator);
    let bxxca = scale_expansion_zeroelim(&bxca, bdx, allocator);
    let byca = scale_expansion_zeroelim(&ca, bdy, allocator);
    let byyca = scale_expansion_zeroelim(&byca, bdy, allocator);
    let bdet = fast_expansion_sum_zeroelim(&bxxca, &byyca, allocator);

    let (adxbdy1, adxbdy0) = two_product(adx, bdy);
    let (bdxady1, bdxady0) = two_product(bdx, ady);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let cxab = scale_expansion_zeroelim(&ab, cdx, allocator);
    let cxxab = scale_expansion_zeroelim(&cxab, cdx, allocator);
    let cyab = scale_expansion_zeroelim(&ab, cdy, allocator);
    let cyyab = scale_expansion_zeroelim(&cyab, cdy, allocator);
    let cdet = fast_expansion_sum_zeroelim(&cxxab, &cyyab, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let mut fin = fast_expansion_sum_zeroelim(&abdet, &cdet, allocator);

    let mut det = estimate(&fin);
    let errbound = B.icc_err_boundb * permanent;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let adxtail = two_diff_tail(pa[0], pd[0], adx);
    let adytail = two_diff_tail(pa[1], pd[1], ady);
    let bdxtail = two_diff_tail(pb[0], pd[0], bdx);
    let bdytail = two_diff_tail(pb[1], pd[1], bdy);
    let cdxtail = two_diff_tail(pc[0], pd[0], cdx);
    let cdytail = two_diff_tail(pc[1], pd[1], cdy);
    if (adxtail == 0.0)
        && (bdxtail == 0.0)
        && (cdxtail == 0.0)
        && (adytail == 0.0)
        && (bdytail == 0.0)
        && (cdytail == 0.0)
    {
        return det;
    }

    let errbound = B.icc_err_boundc * permanent + B.result_err_bound * det.abs();
    det += ((adx * adx + ady * ady)
        * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
        + ((bdx * bdx + bdy * bdy)
            * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
        + ((cdx * cdx + cdy * cdy)
            * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let mut aa = Vec::new_in(allocator);
    let mut bb = Vec::new_in(allocator);
    let mut cc = Vec::new_in(allocator);
    if (bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0) {
        let (adxadx1, adxadx0) = square(adx);
        let (adyady1, adyady0) = square(ady);
        aa.resize(4, 0.0);
        (aa[3], aa[2], aa[1], aa[0]) = two_two_sum(adxadx1, adxadx0, adyady1, adyady0);
    }
    if (cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0) {
        let (bdxbdx1, bdxbdx0) = square(bdx);
        let (bdybdy1, bdybdy0) = square(bdy);
        bb.resize(4, 0.0);
        (bb[3], bb[2], bb[1], bb[0]) = two_two_sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0);
    }
    if (adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0) {
        let (cdxcdx1, cdxcdx0) = square(cdx);
        let (cdycdy1, cdycdy0) = square(cdy);
        cc.resize(4, 0.0);
        (cc[3], cc[2], cc[1], cc[0]) = two_two_sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0);
    }

    fn helper1<A: Allocator + Copy>(
        adxtail: f64,
        bc: &[f64],
        adx: f64,
        cc: &[f64],
        bdy: f64,
        bb: &[f64],
        cdy: f64,
        fin: &mut Vec<f64, A>,
        allocator: A,
    ) -> Vec<f64, A> {
        let axtbc = scale_expansion_zeroelim(&bc, adxtail, allocator);
        let temp16a = scale_expansion_zeroelim(&axtbc, 2.0 * adx, allocator);

        let axtcc = scale_expansion_zeroelim(&cc, adxtail, allocator);
        let temp16b = scale_expansion_zeroelim(&axtcc, bdy, allocator);

        let axtbb = scale_expansion_zeroelim(&bb, adxtail, allocator);
        let temp16c = scale_expansion_zeroelim(&axtbb, -cdy, allocator);

        let temp32a = fast_expansion_sum_zeroelim(&temp16a, &temp16b, allocator);
        let temp48 = fast_expansion_sum_zeroelim(&temp16c, &temp32a, allocator);
        *fin = fast_expansion_sum_zeroelim(&fin, &temp48, allocator);

        axtbc
    }
    let mut axtbc = Vec::new_in(allocator);
    if adxtail != 0.0 {
        axtbc = helper1(adxtail, &bc, adx, &cc, bdy, &bb, cdy, &mut fin, allocator);
    }

    let mut aytbc = Vec::new_in(allocator);
    if adytail != 0.0 {
        aytbc = helper1(adytail, &bc, ady, &bb, cdx, &cc, bdx, &mut fin, allocator);
    }

    let mut bxtca = Vec::new_in(allocator);
    if bdxtail != 0.0 {
        bxtca = helper1(bdxtail, &ca, bdx, &aa, cdy, &cc, ady, &mut fin, allocator);
    }

    let mut bytca = Vec::new_in(allocator);
    if bdytail != 0.0 {
        bytca = helper1(bdytail, &ca, bdy, &cc, adx, &aa, cdx, &mut fin, allocator);
    }

    let mut cxtab = Vec::new_in(allocator);
    if cdxtail != 0.0 {
        cxtab = helper1(cdxtail, &ab, cdx, &bb, ady, &aa, bdy, &mut fin, allocator);
    }

    let mut cytab = Vec::new_in(allocator);
    if cdytail != 0.0 {
        cytab = helper1(cdytail, &ab, cdy, &aa, bdx, &bb, adx, &mut fin, allocator);
    }

    fn helper2<A: Allocator + Copy>(
        adxtail: f64,
        adytail: f64,
        bdxtail: f64,
        bdytail: f64,
        cdxtail: f64,
        cdytail: f64,
        adx: f64,
        ady: f64,
        bdx: f64,
        bdy: f64,
        cdx: f64,
        cdy: f64,
        axtbc: &[f64],
        aytbc: &[f64],
        bb: &[f64],
        cc: &[f64],
        fin: &mut Vec<f64, A>,
        allocator: A,
    ) {
        if (adxtail != 0.0) || (adytail != 0.0) {
            let (bct, bctt) =
                if (bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0) {
                    let (ti1, ti0) = two_product(bdxtail, cdy);
                    let (tj1, tj0) = two_product(bdx, cdytail);
                    let mut u = std::vec::from_elem_in(0.0, 4, allocator);
                    (u[3], u[2], u[1], u[0]) = two_two_sum(ti1, ti0, tj1, tj0);
                    let (ti1, ti0) = two_product(cdxtail, -bdy);
                    let (tj1, tj0) = two_product(cdx, -bdytail);
                    let mut v = std::vec::from_elem_in(0.0, 4, allocator);
                    (v[3], v[2], v[1], v[0]) = two_two_sum(ti1, ti0, tj1, tj0);
                    let bct = fast_expansion_sum_zeroelim(&u, &v, allocator);

                    let (ti1, ti0) = two_product(bdxtail, cdytail);
                    let (tj1, tj0) = two_product(cdxtail, bdytail);
                    let mut bctt = std::vec::from_elem_in(0.0, 4, allocator);
                    (bctt[3], bctt[2], bctt[1], bctt[0]) = two_two_diff(ti1, ti0, tj1, tj0);
                    (bct, bctt)
                } else {
                    ([0.0].to_vec_in(allocator), [0.0].to_vec_in(allocator))
                };

            if adxtail != 0.0 {
                let temp16a = scale_expansion_zeroelim(axtbc, adxtail, allocator);
                let axtbct = scale_expansion_zeroelim(&bct, adxtail, allocator);
                let temp32a = scale_expansion_zeroelim(&axtbct, 2.0 * adx, allocator);
                let temp48 = fast_expansion_sum_zeroelim(&temp16a, &temp32a, allocator);
                *fin = fast_expansion_sum_zeroelim(&fin, &temp48, allocator);
                if bdytail != 0.0 {
                    let temp8 = scale_expansion_zeroelim(cc, adxtail, allocator);
                    let temp16a = scale_expansion_zeroelim(&temp8, bdytail, allocator);
                    *fin = fast_expansion_sum_zeroelim(&fin, &temp16a, allocator);
                }
                if cdytail != 0.0 {
                    let temp8 = scale_expansion_zeroelim(bb, -adxtail, allocator);
                    let temp16a = scale_expansion_zeroelim(&temp8, cdytail, allocator);
                    *fin = fast_expansion_sum_zeroelim(&fin, &temp16a, allocator);
                }

                let temp32a = scale_expansion_zeroelim(&axtbct, adxtail, allocator);
                let axtbctt = scale_expansion_zeroelim(&bctt, adxtail, allocator);
                let temp16a = scale_expansion_zeroelim(&axtbctt, 2.0 * adx, allocator);
                let temp16b = scale_expansion_zeroelim(&axtbctt, adxtail, allocator);
                let temp32b = fast_expansion_sum_zeroelim(&temp16a, &temp16b, allocator);
                let temp64 = fast_expansion_sum_zeroelim(&temp32a, &temp32b, allocator);
                *fin = fast_expansion_sum_zeroelim(&fin, &temp64, allocator);
            }
            if adytail != 0.0 {
                let temp16a = scale_expansion_zeroelim(aytbc, adytail, allocator);
                let aytbct = scale_expansion_zeroelim(&bct, adytail, allocator);
                let temp32a = scale_expansion_zeroelim(&aytbct, 2.0 * ady, allocator);
                let temp48 = fast_expansion_sum_zeroelim(&temp16a, &temp32a, allocator);
                *fin = fast_expansion_sum_zeroelim(&fin, &temp48, allocator);

                let temp32a = scale_expansion_zeroelim(&aytbct, adytail, allocator);
                let aytbctt = scale_expansion_zeroelim(&bctt, adytail, allocator);
                let temp16a = scale_expansion_zeroelim(&aytbctt, 2.0 * ady, allocator);
                let temp16b = scale_expansion_zeroelim(&aytbctt, adytail, allocator);
                let temp32b = fast_expansion_sum_zeroelim(&temp16a, &temp16b, allocator);
                let temp64 = fast_expansion_sum_zeroelim(&temp32a, &temp32b, allocator);
                *fin = fast_expansion_sum_zeroelim(&fin, &temp64, allocator);
            }
        }
    }
    helper2(
        adxtail, adytail, bdxtail, bdytail, cdxtail, cdytail, adx, ady, bdx, bdy, cdx, cdy, &axtbc,
        &aytbc, &bb, &cc, &mut fin, allocator,
    );
    helper2(
        bdxtail, bdytail, cdxtail, cdytail, adxtail, adytail, bdx, bdy, cdx, cdy, adx, ady, &bxtca,
        &bytca, &cc, &aa, &mut fin, allocator,
    );
    helper2(
        cdxtail, cdytail, adxtail, adytail, bdxtail, bdytail, cdx, cdy, adx, ady, bdx, bdy, &cxtab,
        &cytab, &aa, &bb, &mut fin, allocator,
    );
    *fin.last().unwrap()
}

pub fn incircle<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    allocator: A,
) -> f64 {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;
    let errbound = B.icc_err_bounda * permanent;
    if (det > errbound) || (-det > errbound) {
        return det;
    }

    incircle_adapt(pa, pb, pc, pd, permanent, allocator)
}

fn insphere_exact<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    allocator: A,
) -> f64 {
    let (axby1, axby0) = two_product(pa[0], pb[1]);
    let (bxay1, bxay0) = two_product(pb[0], pa[1]);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(axby1, axby0, bxay1, bxay0);

    let (bxcy1, bxcy0) = two_product(pb[0], pc[1]);
    let (cxby1, cxby0) = two_product(pc[0], pb[1]);
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);

    let (cxdy1, cxdy0) = two_product(pc[0], pd[1]);
    let (dxcy1, dxcy0) = two_product(pd[0], pc[1]);
    let mut cd = std::vec::from_elem_in(0.0, 4, allocator);
    (cd[3], cd[2], cd[1], cd[0]) = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);

    let (dxey1, dxey0) = two_product(pd[0], pe[1]);
    let (exdy1, exdy0) = two_product(pe[0], pd[1]);
    let mut de = std::vec::from_elem_in(0.0, 4, allocator);
    (de[3], de[2], de[1], de[0]) = two_two_diff(dxey1, dxey0, exdy1, exdy0);

    let (exay1, exay0) = two_product(pe[0], pa[1]);
    let (axey1, axey0) = two_product(pa[0], pe[1]);
    let mut ea = std::vec::from_elem_in(0.0, 4, allocator);
    (ea[3], ea[2], ea[1], ea[0]) = two_two_diff(exay1, exay0, axey1, axey0);

    let (axcy1, axcy0) = two_product(pa[0], pc[1]);
    let (cxay1, cxay0) = two_product(pc[0], pa[1]);
    let mut ac = std::vec::from_elem_in(0.0, 4, allocator);
    (ac[3], ac[2], ac[1], ac[0]) = two_two_diff(axcy1, axcy0, cxay1, cxay0);

    let (bxdy1, bxdy0) = two_product(pb[0], pd[1]);
    let (dxby1, dxby0) = two_product(pd[0], pb[1]);
    let mut bd = std::vec::from_elem_in(0.0, 4, allocator);
    (bd[3], bd[2], bd[1], bd[0]) = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);

    let (cxey1, cxey0) = two_product(pc[0], pe[1]);
    let (excy1, excy0) = two_product(pe[0], pc[1]);
    let mut ce = std::vec::from_elem_in(0.0, 4, allocator);
    (ce[3], ce[2], ce[1], ce[0]) = two_two_diff(cxey1, cxey0, excy1, excy0);

    let (dxay1, dxay0) = two_product(pd[0], pa[1]);
    let (axdy1, axdy0) = two_product(pa[0], pd[1]);
    let mut da = std::vec::from_elem_in(0.0, 4, allocator);
    (da[3], da[2], da[1], da[0]) = two_two_diff(dxay1, dxay0, axdy1, axdy0);

    let (exby1, exby0) = two_product(pe[0], pb[1]);
    let (bxey1, bxey0) = two_product(pb[0], pe[1]);
    let mut eb = std::vec::from_elem_in(0.0, 4, allocator);
    (eb[3], eb[2], eb[1], eb[0]) = two_two_diff(exby1, exby0, bxey1, bxey0);

    let temp8a = scale_expansion_zeroelim(&bc, pa[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ac, -pb[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ab, pc[2], allocator);
    let abc = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&cd, pb[2], allocator);
    let temp8b = scale_expansion_zeroelim(&bd, -pc[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&bc, pd[2], allocator);
    let bcd = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&de, pc[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ce, -pd[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&cd, pe[2], allocator);
    let cde = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ea, pd[2], allocator);
    let temp8b = scale_expansion_zeroelim(&da, -pe[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&de, pa[2], allocator);
    let dea = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ab, pe[2], allocator);
    let temp8b = scale_expansion_zeroelim(&eb, -pa[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ea, pb[2], allocator);
    let eab = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&bd, pa[2], allocator);
    let temp8b = scale_expansion_zeroelim(&da, pb[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ab, pd[2], allocator);
    let abd = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ce, pb[2], allocator);
    let temp8b = scale_expansion_zeroelim(&eb, pc[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&bc, pe[2], allocator);
    let bce = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&da, pc[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ac, pd[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&cd, pa[2], allocator);
    let cda = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&eb, pd[2], allocator);
    let temp8b = scale_expansion_zeroelim(&bd, pe[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&de, pb[2], allocator);
    let deb = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ac, pe[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ce, pa[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ea, pc[2], allocator);
    let eac = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&cde, &bce, allocator);
    let mut temp48b = fast_expansion_sum_zeroelim(&deb, &bcd, allocator);
    expansion_invert(&mut temp48b);

    let bcde = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let temp192 = scale_expansion_zeroelim(&bcde, pa[0], allocator);
    let det384x = scale_expansion_zeroelim(&temp192, pa[0], allocator);
    let temp192 = scale_expansion_zeroelim(&bcde, pa[1], allocator);
    let det384y = scale_expansion_zeroelim(&temp192, pa[1], allocator);
    let temp192 = scale_expansion_zeroelim(&bcde, pa[2], allocator);
    let det384z = scale_expansion_zeroelim(&temp192, pa[2], allocator);
    let detxy = fast_expansion_sum_zeroelim(&det384x, &det384y, allocator);
    let adet = fast_expansion_sum_zeroelim(&detxy, &det384z, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&dea, &cda, allocator);
    temp48b = fast_expansion_sum_zeroelim(&eac, &cde, allocator);
    expansion_invert(&mut temp48b);

    let cdea = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let temp192 = scale_expansion_zeroelim(&cdea, pb[0], allocator);
    let det384x = scale_expansion_zeroelim(&temp192, pb[0], allocator);
    let temp192 = scale_expansion_zeroelim(&cdea, pb[1], allocator);
    let det384y = scale_expansion_zeroelim(&temp192, pb[1], allocator);
    let temp192 = scale_expansion_zeroelim(&cdea, pb[2], allocator);
    let det384z = scale_expansion_zeroelim(&temp192, pb[2], allocator);
    let detxy = fast_expansion_sum_zeroelim(&det384x, &det384y, allocator);
    let bdet = fast_expansion_sum_zeroelim(&detxy, &det384z, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&eab, &deb, allocator);
    temp48b = fast_expansion_sum_zeroelim(&abd, &dea, allocator);
    expansion_invert(&mut temp48b);

    let deab = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let temp192 = scale_expansion_zeroelim(&deab, pc[0], allocator);
    let det384x = scale_expansion_zeroelim(&temp192, pc[0], allocator);
    let temp192 = scale_expansion_zeroelim(&deab, pc[1], allocator);
    let det384y = scale_expansion_zeroelim(&temp192, pc[1], allocator);
    let temp192 = scale_expansion_zeroelim(&deab, pc[2], allocator);
    let det384z = scale_expansion_zeroelim(&temp192, pc[2], allocator);
    let detxy = fast_expansion_sum_zeroelim(&det384x, &det384y, allocator);
    let cdet = fast_expansion_sum_zeroelim(&detxy, &det384z, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&abc, &eac, allocator);
    temp48b = fast_expansion_sum_zeroelim(&bce, &eab, allocator);
    expansion_invert(&mut temp48b);

    let eabc = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let temp192 = scale_expansion_zeroelim(&eabc, pd[0], allocator);
    let det384x = scale_expansion_zeroelim(&temp192, pd[0], allocator);
    let temp192 = scale_expansion_zeroelim(&eabc, pd[1], allocator);
    let det384y = scale_expansion_zeroelim(&temp192, pd[1], allocator);
    let temp192 = scale_expansion_zeroelim(&eabc, pd[2], allocator);
    let det384z = scale_expansion_zeroelim(&temp192, pd[2], allocator);
    let detxy = fast_expansion_sum_zeroelim(&det384x, &det384y, allocator);
    let ddet = fast_expansion_sum_zeroelim(&detxy, &det384z, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&bcd, &abd, allocator);
    temp48b = fast_expansion_sum_zeroelim(&cda, &abc, allocator);
    expansion_invert(&mut temp48b);

    let abcd = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let temp192 = scale_expansion_zeroelim(&abcd, pe[0], allocator);
    let det384x = scale_expansion_zeroelim(&temp192, pe[0], allocator);
    let temp192 = scale_expansion_zeroelim(&abcd, pe[1], allocator);
    let det384y = scale_expansion_zeroelim(&temp192, pe[1], allocator);
    let temp192 = scale_expansion_zeroelim(&abcd, pe[2], allocator);
    let det384z = scale_expansion_zeroelim(&temp192, pe[2], allocator);
    let detxy = fast_expansion_sum_zeroelim(&det384x, &det384y, allocator);
    let edet = fast_expansion_sum_zeroelim(&detxy, &det384z, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let cddet = fast_expansion_sum_zeroelim(&cdet, &ddet, allocator);
    let cdedet = fast_expansion_sum_zeroelim(&cddet, &edet, allocator);
    let deter = fast_expansion_sum_zeroelim(&abdet, &cdedet, allocator);

    *deter.last().unwrap()
}

fn insphere_adapt<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    permanent: f64,
    allocator: A,
) -> f64 {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];

    let (aexbey1, aexbey0) = two_product(aex, bey);
    let (bexaey1, bexaey0) = two_product(bex, aey);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(aexbey1, aexbey0, bexaey1, bexaey0);

    let (bexcey1, bexcey0) = two_product(bex, cey);
    let (cexbey1, cexbey0) = two_product(cex, bey);
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bexcey1, bexcey0, cexbey1, cexbey0);

    let (cexdey1, cexdey0) = two_product(cex, dey);
    let (dexcey1, dexcey0) = two_product(dex, cey);
    let mut cd = std::vec::from_elem_in(0.0, 4, allocator);
    (cd[3], cd[2], cd[1], cd[0]) = two_two_diff(cexdey1, cexdey0, dexcey1, dexcey0);

    let (dexaey1, dexaey0) = two_product(dex, aey);
    let (aexdey1, aexdey0) = two_product(aex, dey);
    let mut da = std::vec::from_elem_in(0.0, 4, allocator);
    (da[3], da[2], da[1], da[0]) = two_two_diff(dexaey1, dexaey0, aexdey1, aexdey0);

    let (aexcey1, aexcey0) = two_product(aex, cey);
    let (cexaey1, cexaey0) = two_product(cex, aey);
    let mut ac = std::vec::from_elem_in(0.0, 4, allocator);
    (ac[3], ac[2], ac[1], ac[0]) = two_two_diff(aexcey1, aexcey0, cexaey1, cexaey0);

    let (bexdey1, bexdey0) = two_product(bex, dey);
    let (dexbey1, dexbey0) = two_product(dex, bey);
    let mut bd = std::vec::from_elem_in(0.0, 4, allocator);
    (bd[3], bd[2], bd[1], bd[0]) = two_two_diff(bexdey1, bexdey0, dexbey1, dexbey0);

    let temp8a = scale_expansion_zeroelim(&cd, bez, allocator);
    let temp8b = scale_expansion_zeroelim(&bd, -cez, allocator);
    let temp8c = scale_expansion_zeroelim(&bc, dez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, aex, allocator);
    let xdet = scale_expansion_zeroelim(&temp48, -aex, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, aey, allocator);
    let ydet = scale_expansion_zeroelim(&temp48, -aey, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, aez, allocator);
    let zdet = scale_expansion_zeroelim(&temp48, -aez, allocator);
    let xydet = fast_expansion_sum_zeroelim(&xdet, &ydet, allocator);
    let adet = fast_expansion_sum_zeroelim(&xydet, &zdet, allocator);

    let temp8a = scale_expansion_zeroelim(&da, cez, allocator);
    let temp8b = scale_expansion_zeroelim(&ac, dez, allocator);
    let temp8c = scale_expansion_zeroelim(&cd, aez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, bex, allocator);
    let xdet = scale_expansion_zeroelim(&temp48, bex, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, bey, allocator);
    let ydet = scale_expansion_zeroelim(&temp48, bey, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, bez, allocator);
    let zdet = scale_expansion_zeroelim(&temp48, bez, allocator);
    let xydet = fast_expansion_sum_zeroelim(&xdet, &ydet, allocator);
    let bdet = fast_expansion_sum_zeroelim(&xydet, &zdet, allocator);

    let temp8a = scale_expansion_zeroelim(&ab, dez, allocator);
    let temp8b = scale_expansion_zeroelim(&bd, aez, allocator);
    let temp8c = scale_expansion_zeroelim(&da, bez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, cex, allocator);
    let xdet = scale_expansion_zeroelim(&temp48, -cex, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, cey, allocator);
    let ydet = scale_expansion_zeroelim(&temp48, -cey, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, cez, allocator);
    let zdet = scale_expansion_zeroelim(&temp48, -cez, allocator);
    let xydet = fast_expansion_sum_zeroelim(&xdet, &ydet, allocator);
    let cdet = fast_expansion_sum_zeroelim(&xydet, &zdet, allocator);

    let temp8a = scale_expansion_zeroelim(&bc, aez, allocator);
    let temp8b = scale_expansion_zeroelim(&ac, -bez, allocator);
    let temp8c = scale_expansion_zeroelim(&ab, cez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, dex, allocator);
    let xdet = scale_expansion_zeroelim(&temp48, dex, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, dey, allocator);
    let ydet = scale_expansion_zeroelim(&temp48, dey, allocator);
    let temp48 = scale_expansion_zeroelim(&temp24, dez, allocator);
    let zdet = scale_expansion_zeroelim(&temp48, dez, allocator);
    let xydet = fast_expansion_sum_zeroelim(&xdet, &ydet, allocator);
    let ddet = fast_expansion_sum_zeroelim(&xydet, &zdet, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let cddet = fast_expansion_sum_zeroelim(&cdet, &ddet, allocator);
    let fin1 = fast_expansion_sum_zeroelim(&abdet, &cddet, allocator);

    let mut det = estimate(&fin1);
    let errbound = B.isp_err_boundb * permanent;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let aextail = two_diff_tail(pa[0], pe[0], aex);
    let aeytail = two_diff_tail(pa[1], pe[1], aey);
    let aeztail = two_diff_tail(pa[2], pe[2], aez);
    let bextail = two_diff_tail(pb[0], pe[0], bex);
    let beytail = two_diff_tail(pb[1], pe[1], bey);
    let beztail = two_diff_tail(pb[2], pe[2], bez);
    let cextail = two_diff_tail(pc[0], pe[0], cex);
    let ceytail = two_diff_tail(pc[1], pe[1], cey);
    let ceztail = two_diff_tail(pc[2], pe[2], cez);
    let dextail = two_diff_tail(pd[0], pe[0], dex);
    let deytail = two_diff_tail(pd[1], pe[1], dey);
    let deztail = two_diff_tail(pd[2], pe[2], dez);
    if (aextail == 0.0)
        && (aeytail == 0.0)
        && (aeztail == 0.0)
        && (bextail == 0.0)
        && (beytail == 0.0)
        && (beztail == 0.0)
        && (cextail == 0.0)
        && (ceytail == 0.0)
        && (ceztail == 0.0)
        && (dextail == 0.0)
        && (deytail == 0.0)
        && (deztail == 0.0)
    {
        return det;
    }

    let errbound = B.isp_err_boundc * permanent + B.result_err_bound * det.abs();
    let abeps = (aex * beytail + bey * aextail) - (aey * bextail + bex * aeytail);
    let bceps = (bex * ceytail + cey * bextail) - (bey * cextail + cex * beytail);
    let cdeps = (cex * deytail + dey * cextail) - (cey * dextail + dex * ceytail);
    let daeps = (dex * aeytail + aey * dextail) - (dey * aextail + aex * deytail);
    let aceps = (aex * ceytail + cey * aextail) - (aey * cextail + cex * aeytail);
    let bdeps = (bex * deytail + dey * bextail) - (bey * dextail + dex * beytail);
    det += (((bex * bex + bey * bey + bez * bez)
        * ((cez * daeps + dez * aceps + aez * cdeps)
            + (ceztail * da[3] + deztail * ac[3] + aeztail * cd[3]))
        + (dex * dex + dey * dey + dez * dez)
            * ((aez * bceps - bez * aceps + cez * abeps)
                + (aeztail * bc[3] - beztail * ac[3] + ceztail * ab[3])))
        - ((aex * aex + aey * aey + aez * aez)
            * ((bez * cdeps - cez * bdeps + dez * bceps)
                + (beztail * cd[3] - ceztail * bd[3] + deztail * bc[3]))
            + (cex * cex + cey * cey + cez * cez)
                * ((dez * abeps + aez * bdeps + bez * daeps)
                    + (deztail * ab[3] + aeztail * bd[3] + beztail * da[3]))))
        + 2.0
            * (((bex * bextail + bey * beytail + bez * beztail)
                * (cez * da[3] + dez * ac[3] + aez * cd[3])
                + (dex * dextail + dey * deytail + dez * deztail)
                    * (aez * bc[3] - bez * ac[3] + cez * ab[3]))
                - ((aex * aextail + aey * aeytail + aez * aeztail)
                    * (bez * cd[3] - cez * bd[3] + dez * bc[3])
                    + (cex * cextail + cey * ceytail + cez * ceztail)
                        * (dez * ab[3] + aez * bd[3] + bez * da[3])));
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    insphere_exact(pa, pb, pc, pd, pe, allocator)
}

pub fn insphere<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    allocator: A,
) -> f64 {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];

    let aexbey = aex * bey;
    let bexaey = bex * aey;
    let ab = aexbey - bexaey;
    let bexcey = bex * cey;
    let cexbey = cex * bey;
    let bc = bexcey - cexbey;
    let cexdey = cex * dey;
    let dexcey = dex * cey;
    let cd = cexdey - dexcey;
    let dexaey = dex * aey;
    let aexdey = aex * dey;
    let da = dexaey - aexdey;

    let aexcey = aex * cey;
    let cexaey = cex * aey;
    let ac = aexcey - cexaey;
    let bexdey = bex * dey;
    let dexbey = dex * bey;
    let bd = bexdey - dexbey;

    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;

    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;

    let det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    let max_var = aex
        .abs()
        .max(bex.abs())
        .max(cex.abs())
        .max(dex.abs())
        .max(aey.abs())
        .max(bey.abs())
        .max(cey.abs())
        .max(dey.abs())
        .max(aez.abs())
        .max(bez.abs())
        .max(cez.abs())
        .max(dez.abs());
    let mut epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= epsilon;
    epsilon *= max_var;
    epsilon *= 1.145750161413163e-13;
    if (det > epsilon) || (-det > epsilon) {
        return det;
    }

    let aezplus = aez.abs();
    let bezplus = bez.abs();
    let cezplus = cez.abs();
    let dezplus = dez.abs();
    let aexbeyplus = aexbey.abs();
    let bexaeyplus = bexaey.abs();
    let bexceyplus = bexcey.abs();
    let cexbeyplus = cexbey.abs();
    let cexdeyplus = cexdey.abs();
    let dexceyplus = dexcey.abs();
    let dexaeyplus = dexaey.abs();
    let aexdeyplus = aexdey.abs();
    let aexceyplus = aexcey.abs();
    let cexaeyplus = cexaey.abs();
    let bexdeyplus = bexdey.abs();
    let dexbeyplus = dexbey.abs();
    let permanent = ((cexdeyplus + dexceyplus) * bezplus
        + (dexbeyplus + bexdeyplus) * cezplus
        + (bexceyplus + cexbeyplus) * dezplus)
        * alift
        + ((dexaeyplus + aexdeyplus) * cezplus
            + (aexceyplus + cexaeyplus) * dezplus
            + (cexdeyplus + dexceyplus) * aezplus)
            * blift
        + ((aexbeyplus + bexaeyplus) * dezplus
            + (bexdeyplus + dexbeyplus) * aezplus
            + (dexaeyplus + aexdeyplus) * bezplus)
            * clift
        + ((bexceyplus + cexbeyplus) * aezplus
            + (cexaeyplus + aexceyplus) * bezplus
            + (aexbeyplus + bexaeyplus) * cezplus)
            * dlift;

    insphere_adapt(pa, pb, pc, pd, pe, permanent, allocator)
}

fn orient4d_exact<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    aheight: f64,
    bheight: f64,
    cheight: f64,
    dheight: f64,
    eheight: f64,
    allocator: A,
) -> f64 {
    let (axby1, axby0) = two_product(pa[0], pb[1]);
    let (bxay1, bxay0) = two_product(pb[0], pa[1]);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(axby1, axby0, bxay1, bxay0);

    let (bxcy1, bxcy0) = two_product(pb[0], pc[1]);
    let (cxby1, cxby0) = two_product(pc[0], pb[1]);
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);

    let (cxdy1, cxdy0) = two_product(pc[0], pd[1]);
    let (dxcy1, dxcy0) = two_product(pd[0], pc[1]);
    let mut cd = std::vec::from_elem_in(0.0, 4, allocator);
    (cd[3], cd[2], cd[1], cd[0]) = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);

    let (dxey1, dxey0) = two_product(pd[0], pe[1]);
    let (exdy1, exdy0) = two_product(pe[0], pd[1]);
    let mut de = std::vec::from_elem_in(0.0, 4, allocator);
    (de[3], de[2], de[1], de[0]) = two_two_diff(dxey1, dxey0, exdy1, exdy0);

    let (exay1, exay0) = two_product(pe[0], pa[1]);
    let (axey1, axey0) = two_product(pa[0], pe[1]);
    let mut ea = std::vec::from_elem_in(0.0, 4, allocator);
    (ea[3], ea[2], ea[1], ea[0]) = two_two_diff(exay1, exay0, axey1, axey0);

    let (axcy1, axcy0) = two_product(pa[0], pc[1]);
    let (cxay1, cxay0) = two_product(pc[0], pa[1]);
    let mut ac = std::vec::from_elem_in(0.0, 4, allocator);
    (ac[3], ac[2], ac[1], ac[0]) = two_two_diff(axcy1, axcy0, cxay1, cxay0);

    let (bxdy1, bxdy0) = two_product(pb[0], pd[1]);
    let (dxby1, dxby0) = two_product(pd[0], pb[1]);
    let mut bd = std::vec::from_elem_in(0.0, 4, allocator);
    (bd[3], bd[2], bd[1], bd[0]) = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);

    let (cxey1, cxey0) = two_product(pc[0], pe[1]);
    let (excy1, excy0) = two_product(pe[0], pc[1]);
    let mut ce = std::vec::from_elem_in(0.0, 4, allocator);
    (ce[3], ce[2], ce[1], ce[0]) = two_two_diff(cxey1, cxey0, excy1, excy0);

    let (dxay1, dxay0) = two_product(pd[0], pa[1]);
    let (axdy1, axdy0) = two_product(pa[0], pd[1]);
    let mut da = std::vec::from_elem_in(0.0, 4, allocator);
    (da[3], da[2], da[1], da[0]) = two_two_diff(dxay1, dxay0, axdy1, axdy0);

    let (exby1, exby0) = two_product(pe[0], pb[1]);
    let (bxey1, bxey0) = two_product(pb[0], pe[1]);
    let mut eb = std::vec::from_elem_in(0.0, 4, allocator);
    (eb[3], eb[2], eb[1], eb[0]) = two_two_diff(exby1, exby0, bxey1, bxey0);

    let temp8a = scale_expansion_zeroelim(&bc, pa[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ac, -pb[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ab, pc[2], allocator);
    let abc = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&cd, pb[2], allocator);
    let temp8b = scale_expansion_zeroelim(&bd, -pc[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&bc, pd[2], allocator);
    let bcd = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&de, pc[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ce, -pd[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&cd, pe[2], allocator);
    let cde = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ea, pd[2], allocator);
    let temp8b = scale_expansion_zeroelim(&da, -pe[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&de, pa[2], allocator);
    let dea = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ab, pe[2], allocator);
    let temp8b = scale_expansion_zeroelim(&eb, -pa[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ea, pb[2], allocator);
    let eab = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&bd, pa[2], allocator);
    let temp8b = scale_expansion_zeroelim(&da, pb[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ab, pd[2], allocator);
    let abd = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ce, pb[2], allocator);
    let temp8b = scale_expansion_zeroelim(&eb, pc[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&bc, pe[2], allocator);
    let bce = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&da, pc[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ac, pd[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&cd, pa[2], allocator);
    let cda = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&eb, pd[2], allocator);
    let temp8b = scale_expansion_zeroelim(&bd, pe[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&de, pb[2], allocator);
    let deb = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp8a = scale_expansion_zeroelim(&ac, pe[2], allocator);
    let temp8b = scale_expansion_zeroelim(&ce, pa[2], allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp8a = scale_expansion_zeroelim(&ea, pc[2], allocator);
    let eac = fast_expansion_sum_zeroelim(&temp8a, &temp16, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&cde, &bce, allocator);
    let mut temp48b = fast_expansion_sum_zeroelim(&deb, &bcd, allocator);
    expansion_invert(&mut temp48b);

    let bcde = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let adet = scale_expansion_zeroelim(&bcde, aheight, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&dea, &cda, allocator);
    temp48b = fast_expansion_sum_zeroelim(&eac, &cde, allocator);
    expansion_invert(&mut temp48b);
    let cdea = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let bdet = scale_expansion_zeroelim(&cdea, bheight, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&eab, &deb, allocator);
    temp48b = fast_expansion_sum_zeroelim(&abd, &dea, allocator);
    expansion_invert(&mut temp48b);

    let deab = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let cdet = scale_expansion_zeroelim(&deab, cheight, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&abc, &eac, allocator);
    temp48b = fast_expansion_sum_zeroelim(&bce, &eab, allocator);
    expansion_invert(&mut temp48b);

    let eabc = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let ddet = scale_expansion_zeroelim(&eabc, dheight, allocator);

    let temp48a = fast_expansion_sum_zeroelim(&bcd, &abd, allocator);
    temp48b = fast_expansion_sum_zeroelim(&cda, &abc, allocator);
    expansion_invert(&mut temp48b);

    let abcd = fast_expansion_sum_zeroelim(&temp48a, &temp48b, allocator);
    let edet = scale_expansion_zeroelim(&abcd, eheight, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let cddet = fast_expansion_sum_zeroelim(&cdet, &ddet, allocator);
    let cdedet = fast_expansion_sum_zeroelim(&cddet, &edet, allocator);
    let deter = fast_expansion_sum_zeroelim(&abdet, &cdedet, allocator);

    *deter.last().unwrap()
}

fn orient4d_adapt<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    aheight: f64,
    bheight: f64,
    cheight: f64,
    dheight: f64,
    eheight: f64,
    permanent: f64,
    allocator: A,
) -> f64 {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];
    let aeheight = aheight - eheight;
    let beheight = bheight - eheight;
    let ceheight = cheight - eheight;
    let deheight = dheight - eheight;

    let (aexbey1, aexbey0) = two_product(aex, bey);
    let (bexaey1, bexaey0) = two_product(bex, aey);
    let mut ab = std::vec::from_elem_in(0.0, 4, allocator);
    (ab[3], ab[2], ab[1], ab[0]) = two_two_diff(aexbey1, aexbey0, bexaey1, bexaey0);

    let (bexcey1, bexcey0) = two_product(bex, cey);
    let (cexbey1, cexbey0) = two_product(cex, bey);
    let mut bc = std::vec::from_elem_in(0.0, 4, allocator);
    (bc[3], bc[2], bc[1], bc[0]) = two_two_diff(bexcey1, bexcey0, cexbey1, cexbey0);

    let (cexdey1, cexdey0) = two_product(cex, dey);
    let (dexcey1, dexcey0) = two_product(dex, cey);
    let mut cd = std::vec::from_elem_in(0.0, 4, allocator);
    (cd[3], cd[2], cd[1], cd[0]) = two_two_diff(cexdey1, cexdey0, dexcey1, dexcey0);

    let (dexaey1, dexaey0) = two_product(dex, aey);
    let (aexdey1, aexdey0) = two_product(aex, dey);
    let mut da = std::vec::from_elem_in(0.0, 4, allocator);
    (da[3], da[2], da[1], da[0]) = two_two_diff(dexaey1, dexaey0, aexdey1, aexdey0);

    let (aexcey1, aexcey0) = two_product(aex, cey);
    let (cexaey1, cexaey0) = two_product(cex, aey);
    let mut ac = std::vec::from_elem_in(0.0, 4, allocator);
    (ac[3], ac[2], ac[1], ac[0]) = two_two_diff(aexcey1, aexcey0, cexaey1, cexaey0);

    let (bexdey1, bexdey0) = two_product(bex, dey);
    let (dexbey1, dexbey0) = two_product(dex, bey);
    let mut bd = std::vec::from_elem_in(0.0, 4, allocator);
    (bd[3], bd[2], bd[1], bd[0]) = two_two_diff(bexdey1, bexdey0, dexbey1, dexbey0);

    let temp8a = scale_expansion_zeroelim(&cd, bez, allocator);
    let temp8b = scale_expansion_zeroelim(&bd, -cez, allocator);
    let temp8c = scale_expansion_zeroelim(&bc, dez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let adet = scale_expansion_zeroelim(&temp24, -aeheight, allocator);

    let temp8a = scale_expansion_zeroelim(&da, cez, allocator);
    let temp8b = scale_expansion_zeroelim(&ac, dez, allocator);
    let temp8c = scale_expansion_zeroelim(&cd, aez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let bdet = scale_expansion_zeroelim(&temp24, beheight, allocator);

    let temp8a = scale_expansion_zeroelim(&ab, dez, allocator);
    let temp8b = scale_expansion_zeroelim(&bd, aez, allocator);
    let temp8c = scale_expansion_zeroelim(&da, bez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let cdet = scale_expansion_zeroelim(&temp24, -ceheight, allocator);

    let temp8a = scale_expansion_zeroelim(&bc, aez, allocator);
    let temp8b = scale_expansion_zeroelim(&ac, -bez, allocator);
    let temp8c = scale_expansion_zeroelim(&ab, cez, allocator);
    let temp16 = fast_expansion_sum_zeroelim(&temp8a, &temp8b, allocator);
    let temp24 = fast_expansion_sum_zeroelim(&temp8c, &temp16, allocator);
    let ddet = scale_expansion_zeroelim(&temp24, deheight, allocator);

    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet, allocator);
    let cddet = fast_expansion_sum_zeroelim(&cdet, &ddet, allocator);
    let fin1 = fast_expansion_sum_zeroelim(&abdet, &cddet, allocator);

    let mut det = estimate(&fin1);
    let errbound = B.isp_err_boundb * permanent;
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    let aextail = two_diff_tail(pa[0], pe[0], aex);
    let aeytail = two_diff_tail(pa[1], pe[1], aey);
    let aeztail = two_diff_tail(pa[2], pe[2], aez);
    let aeheighttail = two_diff_tail(aheight, eheight, aeheight);
    let bextail = two_diff_tail(pb[0], pe[0], bex);
    let beytail = two_diff_tail(pb[1], pe[1], bey);
    let beztail = two_diff_tail(pb[2], pe[2], bez);
    let beheighttail = two_diff_tail(bheight, eheight, beheight);
    let cextail = two_diff_tail(pc[0], pe[0], cex);
    let ceytail = two_diff_tail(pc[1], pe[1], cey);
    let ceztail = two_diff_tail(pc[2], pe[2], cez);
    let ceheighttail = two_diff_tail(cheight, eheight, ceheight);
    let dextail = two_diff_tail(pd[0], pe[0], dex);
    let deytail = two_diff_tail(pd[1], pe[1], dey);
    let deztail = two_diff_tail(pd[2], pe[2], dez);
    let deheighttail = two_diff_tail(dheight, eheight, deheight);
    if (aextail == 0.0)
        && (aeytail == 0.0)
        && (aeztail == 0.0)
        && (bextail == 0.0)
        && (beytail == 0.0)
        && (beztail == 0.0)
        && (cextail == 0.0)
        && (ceytail == 0.0)
        && (ceztail == 0.0)
        && (dextail == 0.0)
        && (deytail == 0.0)
        && (deztail == 0.0)
        && (aeheighttail == 0.0)
        && (beheighttail == 0.0)
        && (ceheighttail == 0.0)
        && (deheighttail == 0.0)
    {
        return det;
    }

    let errbound = B.isp_err_boundc * permanent + B.result_err_bound * det.abs();
    let abeps = (aex * beytail + bey * aextail) - (aey * bextail + bex * aeytail);
    let bceps = (bex * ceytail + cey * bextail) - (bey * cextail + cex * beytail);
    let cdeps = (cex * deytail + dey * cextail) - (cey * dextail + dex * ceytail);
    let daeps = (dex * aeytail + aey * dextail) - (dey * aextail + aex * deytail);
    let aceps = (aex * ceytail + cey * aextail) - (aey * cextail + cex * aeytail);
    let bdeps = (bex * deytail + dey * bextail) - (bey * dextail + dex * beytail);
    det += ((beheight
        * ((cez * daeps + dez * aceps + aez * cdeps)
            + (ceztail * da[3] + deztail * ac[3] + aeztail * cd[3]))
        + deheight
            * ((aez * bceps - bez * aceps + cez * abeps)
                + (aeztail * bc[3] - beztail * ac[3] + ceztail * ab[3])))
        - (aeheight
            * ((bez * cdeps - cez * bdeps + dez * bceps)
                + (beztail * cd[3] - ceztail * bd[3] + deztail * bc[3]))
            + ceheight
                * ((dez * abeps + aez * bdeps + bez * daeps)
                    + (deztail * ab[3] + aeztail * bd[3] + beztail * da[3]))))
        + ((beheighttail * (cez * da[3] + dez * ac[3] + aez * cd[3])
            + deheighttail * (aez * bc[3] - bez * ac[3] + cez * ab[3]))
            - (aeheighttail * (bez * cd[3] - cez * bd[3] + dez * bc[3])
                + ceheighttail * (dez * ab[3] + aez * bd[3] + bez * da[3])));
    if (det >= errbound) || (-det >= errbound) {
        return det;
    }

    orient4d_exact(
        pa, pb, pc, pd, pe, aheight, bheight, cheight, dheight, eheight, allocator,
    )
}

pub fn orient4d<A: Allocator + Copy>(
    pa: &[f64],
    pb: &[f64],
    pc: &[f64],
    pd: &[f64],
    pe: &[f64],
    allocator: A,
    aheight: f64,
    bheight: f64,
    cheight: f64,
    dheight: f64,
    eheight: f64,
) -> f64 {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];
    let aeheight = aheight - eheight;
    let beheight = bheight - eheight;
    let ceheight = cheight - eheight;
    let deheight = dheight - eheight;

    let aexbey = aex * bey;
    let bexaey = bex * aey;
    let ab = aexbey - bexaey;
    let bexcey = bex * cey;
    let cexbey = cex * bey;
    let bc = bexcey - cexbey;
    let cexdey = cex * dey;
    let dexcey = dex * cey;
    let cd = cexdey - dexcey;
    let dexaey = dex * aey;
    let aexdey = aex * dey;
    let da = dexaey - aexdey;

    let aexcey = aex * cey;
    let cexaey = cex * aey;
    let ac = aexcey - cexaey;
    let bexdey = bex * dey;
    let dexbey = dex * bey;
    let bd = bexdey - dexbey;

    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;

    let det = (deheight * abc - ceheight * dab) + (beheight * cda - aeheight * bcd);

    let aezplus = aez.abs();
    let bezplus = bez.abs();
    let cezplus = cez.abs();
    let dezplus = dez.abs();
    let aexbeyplus = aexbey.abs();
    let bexaeyplus = bexaey.abs();
    let bexceyplus = bexcey.abs();
    let cexbeyplus = cexbey.abs();
    let cexdeyplus = cexdey.abs();
    let dexceyplus = dexcey.abs();
    let dexaeyplus = dexaey.abs();
    let aexdeyplus = aexdey.abs();
    let aexceyplus = aexcey.abs();
    let cexaeyplus = cexaey.abs();
    let bexdeyplus = bexdey.abs();
    let dexbeyplus = dexbey.abs();
    let permanent = ((cexdeyplus + dexceyplus) * bezplus
        + (dexbeyplus + bexdeyplus) * cezplus
        + (bexceyplus + cexbeyplus) * dezplus)
        * aeheight.abs()
        + ((dexaeyplus + aexdeyplus) * cezplus
            + (aexceyplus + cexaeyplus) * dezplus
            + (cexdeyplus + dexceyplus) * aezplus)
            * beheight.abs()
        + ((aexbeyplus + bexaeyplus) * dezplus
            + (bexdeyplus + dexbeyplus) * aezplus
            + (dexaeyplus + aexdeyplus) * bezplus)
            * ceheight.abs()
        + ((bexceyplus + cexbeyplus) * aezplus
            + (cexaeyplus + aexceyplus) * bezplus
            + (aexbeyplus + bexaeyplus) * cezplus)
            * deheight.abs();
    let errbound = B.isp_err_bounda * permanent;
    if (det > errbound) || (-det > errbound) {
        return det;
    }

    return orient4d_adapt(
        pa, pb, pc, pd, pe, aheight, bheight, cheight, dheight, eheight, permanent, allocator,
    );
}

#[test]
fn test_two_arr_mul() {
    let a = [100.0, 0.0000000001];
    let b = a.clone();
    let pro = mul_expansion_zeroelim(&a, &b, std::alloc::Global);
    assert_eq!(
        estimate(&pro),
        10000.000000020002 - 8.0374068511808128e-13 + 3.1019272970715786e-25
    );
}
