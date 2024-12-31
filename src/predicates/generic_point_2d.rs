use std::{
    alloc::Allocator,
    cell::{Ref, RefCell},
};

use super::{
    abs_max, dummy_abs_max, estimate, get_exponent, ExpansionNum, GenericNum, IntervalNumber,
};

#[derive(Clone, PartialEq, Debug)]
pub struct ExplicitPoint2D<'a> {
    pub data: &'a [f64],
}

impl<'a> From<&'a [f64]> for ExplicitPoint2D<'a> {
    fn from(data: &'a [f64]) -> Self {
        ExplicitPoint2D { data }
    }
}

#[derive(Clone)]
pub struct Implicit2DCache<T> {
    pub x: T,
    pub y: T,
    pub d: T,
}

#[inline(always)]
fn copy_exact_cache<A: Allocator + Copy>(
    cache: &Implicit2DCache<ExpansionNum>,
    allocator: A,
) -> Implicit2DCache<ExpansionNum<A>> {
    Implicit2DCache {
        x: cache.x.to_vec_in(allocator).into(),
        y: cache.y.to_vec_in(allocator).into(),
        d: cache.d.to_vec_in(allocator).into(),
    }
}

pub struct ImplicitPointSSI<'a> {
    pub a: &'a [f64],
    b: &'a [f64],
    p: &'a [f64],
    q: &'a [f64],

    ss_filter: RefCell<Option<(Implicit2DCache<f64>, f64)>>,
    d_filter: RefCell<Option<Implicit2DCache<IntervalNumber>>>,
    exact: RefCell<Option<Implicit2DCache<ExpansionNum>>>,
}

impl<'a> ImplicitPointSSI<'a> {
    pub fn new(a: &'a [f64], b: &'a [f64], p: &'a [f64], q: &'a [f64]) -> Self {
        ImplicitPointSSI {
            a,
            b,
            p,
            q,
            ss_filter: RefCell::new(None),
            d_filter: RefCell::new(None),
            exact: RefCell::new(None),
        }
    }

    pub fn ss_filter(&self) -> Option<&(Implicit2DCache<f64>, f64)> {
        if self.ss_filter.borrow().is_some() {
            let filter = self.ss_filter.borrow();
            if filter.as_ref().unwrap().1 == 0.0 {
                return None;
            } else {
                return Ref::leak(filter).as_ref();
            }
        } else {
            let (mut filter, max_var) = ssi_lambda::<true, _, _>(
                self.a[0], self.a[1], self.b[0], self.b[1], self.p[0], self.p[1], self.q[0],
                self.q[1], abs_max,
            );
            let max_var = max_var.unwrap();
            let mut lambda_d_eps = max_var;
            lambda_d_eps *= lambda_d_eps;
            lambda_d_eps *= 8.881784197001252e-16;
            if filter.d > lambda_d_eps || filter.d < -lambda_d_eps {
                if filter.d < 0.0 {
                    filter.x = -filter.x;
                    filter.y = -filter.y;
                    filter.d = -filter.d;
                }
                self.ss_filter.replace(Some((filter, max_var)));
                return Ref::leak(self.ss_filter.borrow()).as_ref();
            } else {
                self.ss_filter.replace(Some((
                    Implicit2DCache {
                        x: 0.0,
                        y: 0.0,
                        d: 0.0,
                    },
                    0.0,
                )));
                return None;
            }
        }
    }

    pub fn d_filter(&self) -> Option<&Implicit2DCache<IntervalNumber>> {
        if self.d_filter.borrow().is_some() {
            let filter = self.d_filter.borrow();
            if filter.as_ref().unwrap().d.not_zero() {
                return Ref::leak(filter).as_ref();
            } else {
                return None;
            }
        } else {
            let (mut filter, _) = ssi_lambda::<false, IntervalNumber, _>(
                self.a[0].into(),
                self.a[1].into(),
                self.b[0].into(),
                self.b[1].into(),
                self.p[0].into(),
                self.p[1].into(),
                self.q[0].into(),
                self.q[1].into(),
                dummy_abs_max,
            );
            if filter.d.negative() {
                filter.x.neg();
                filter.y.neg();
                filter.d.neg();
            }
            self.d_filter.replace(Some(filter));
            if self.d_filter.borrow().as_ref().unwrap().d.not_zero() {
                return Ref::leak(self.d_filter.borrow()).as_ref();
            } else {
                return None;
            }
        }
    }

    pub fn exact<A: Allocator + Copy>(
        &self,
        allocator: A,
    ) -> Option<Implicit2DCache<ExpansionNum<A>>> {
        if self.exact.borrow().is_some() {
            let exact = self.exact.borrow();
            let exact = exact.as_ref().unwrap();
            if exact.d.not_zero() {
                return Some(copy_exact_cache(exact, allocator));
            } else {
                return None;
            }
        } else {
            let (mut exact, _) = ssi_lambda::<false, ExpansionNum, _>(
                vec![self.a[0]].into(),
                vec![self.a[1]].into(),
                vec![self.b[0]].into(),
                vec![self.b[1]].into(),
                vec![self.p[0]].into(),
                vec![self.p[1]].into(),
                vec![self.q[0]].into(),
                vec![self.q[1]].into(),
                dummy_abs_max,
            );
            if exact.d.negative() {
                exact.x.neg();
                exact.y.neg();
                exact.d.neg();
            }
            normalize_lambda2d(&mut exact.x, &mut exact.y, &mut exact.d);

            self.exact.replace(Some(exact));
            if self.exact.borrow().as_ref().unwrap().d.not_zero() {
                return Some(copy_exact_cache(
                    self.exact.borrow().as_ref().unwrap(),
                    allocator,
                ));
            } else {
                return None;
            }
        }
    }
}

fn ssi_lambda<const NEED_MAX: bool, T: GenericNum, F>(
    ax: T,
    ay: T,
    bx: T,
    by: T,
    px: T,
    py: T,
    qx: T,
    qy: T,
    abs_max: F,
) -> (Implicit2DCache<T>, Option<T>)
where
    F: FnOnce(&[T]) -> Option<T>,
{
    let bax = &ax - &bx;
    let bay = &ay - &by;
    let pqx = qx - &px;
    let pqy = qy - &py;
    let pax = ax - &px;
    let pay = ay - py;
    let d01 = &bax * &pqy;
    let d10 = &bay * &pqx;
    let d = d01 - d10;
    let n01 = pax * &pqy;
    let n10 = pay * &pqx;
    let n = n10 - n01;
    let x = &n * &bax;
    let y = n * &bay;
    if NEED_MAX {
        (Implicit2DCache { x, y, d }, abs_max(&[bax, bay, pqx, pqy]))
    } else {
        (Implicit2DCache { x, y, d }, None)
    }
}

fn normalize_lambda2d(x: &mut [f64], y: &mut [f64], d: &mut [f64]) {
    let data = [x, y, d];
    let max_val = data
        .iter()
        .map(|arr| estimate(arr))
        .max_by(|x, y| x.abs().total_cmp(&y.abs()))
        .unwrap();
    let e = get_exponent(max_val);
    if e != 0 {
        let s = 2.0f64.powi(-e);
        for arr in data {
            for val in arr {
                *val *= s;
            }
        }
    }
}

pub enum Point2D<'a> {
    E(&'a [f64]),
    I(ImplicitPointSSI<'a>),
}
