use std::{
    cell::{Ref, RefCell},
    ops::Deref,
};

use bumpalo::{collections::Vec, vec, Bump};

use super::{
    abs_max, dummy_abs_max, estimate, get_exponent, ExpansionNum, GenericNum, IntervalNumber,
};

#[derive(Clone, PartialEq)]
pub struct ExplicitPoint3D {
    pub data: [f64; 3],
}

impl Deref for ExplicitPoint3D {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl From<&[f64]> for ExplicitPoint3D {
    fn from(d: &[f64]) -> Self {
        Self {
            data: [d[0], d[1], d[2]],
        }
    }
}

#[inline(always)]
fn normalize_lambda3d(x: &mut [f64], y: &mut [f64], z: &mut [f64], d: &mut [f64]) {
    let data = [x, y, z, d];
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

#[derive(Clone)]
pub struct Implicit3DCache<T> {
    pub x: T,
    pub y: T,
    pub z: T,
    pub d: T,
}

pub trait ImplicitPoint3D<'b> {
    fn static_filter(&self) -> Option<&(Implicit3DCache<f64>, f64)>;
    fn dynamic_filter(&self) -> Option<&Implicit3DCache<IntervalNumber>>;
    fn exact(&self) -> Option<&Implicit3DCache<ExpansionNum<'b>>>;
}

/// A point in 3D space representing the intersection of a line and a plane.
#[derive(Clone)]
pub struct ImplicitPointLPI<'b> {
    /// line
    pub p: ExplicitPoint3D,
    pub q: ExplicitPoint3D,
    /// plane
    pub r: ExplicitPoint3D,
    pub s: ExplicitPoint3D,
    pub t: ExplicitPoint3D,

    bump: &'b Bump,

    ss_filter: RefCell<Option<(Implicit3DCache<f64>, f64)>>,
    d_filter: RefCell<Option<Implicit3DCache<IntervalNumber>>>,
    exact: RefCell<Option<Implicit3DCache<ExpansionNum<'b>>>>,
}

impl<'b> ImplicitPointLPI<'b> {
    pub fn new(
        p: ExplicitPoint3D,
        q: ExplicitPoint3D,
        r: ExplicitPoint3D,
        s: ExplicitPoint3D,
        t: ExplicitPoint3D,
        bump: &'b Bump,
    ) -> Self {
        Self {
            p,
            q,
            r,
            s,
            t,
            ss_filter: RefCell::new(None),
            d_filter: RefCell::new(None),
            exact: RefCell::new(None),
            bump,
        }
    }
}

fn lpi_lambda<'b, const NEED_MAX: bool, T: GenericNum + 'b, F>(
    px: T,
    py: T,
    pz: T,
    qx: T,
    qy: T,
    qz: T,
    rx: T,
    ry: T,
    rz: T,
    sx: T,
    sy: T,
    sz: T,
    tx: T,
    ty: T,
    tz: T,
    abs_max: F,
    bump: &'b Bump,
) -> (Implicit3DCache<T>, Option<T>)
where
    F: FnOnce(Vec<'b, T>) -> Option<T>,
{
    let a11 = &px - &qx;
    let a12 = &py - &qy;
    let a13 = &pz - &qz;
    let a21 = sx - &rx;
    let a22 = sy - &ry;
    let a23 = sz - &rz;
    let a31 = tx - &rx;
    let a32 = ty - &ry;
    let a33 = tz - &rz;
    let tv1 = &a22 * &a33;
    let tv2 = &a23 * &a32;
    let a2233 = tv1 - tv2;
    let tv3 = &a21 * &a33;
    let tv4 = &a23 * &a31;
    let a2133 = tv3 - tv4;
    let tv5 = &a21 * &a32;
    let tv6 = &a22 * &a31;
    let a2132 = tv5 - tv6;
    let tv7 = &a11 * &a2233;
    let tv8 = &a12 * &a2133;
    let tv9 = &a13 * &a2132;
    let tt1 = tv7 - tv8;

    let d = tt1 + tv9;

    let px_rx = &px - &rx;
    let py_ry = &py - &ry;
    let pz_rz = &pz - &rz;
    let tt2 = &py_ry * &a2133;
    let tt3 = &px_rx * &a2233;
    let tt4 = &pz_rz * &a2132;
    let tt5 = tt3 + tt4;
    let n = tt5 - tt2;
    let ax = &a11 * &n;
    let ay = &a12 * &n;
    let az = &a13 * &n;
    let dpx = &d * &px;
    let dpy = &d * &py;
    let dpz = &d * &pz;
    let x = dpx - ax;
    let y = dpy - ay;
    let z = dpz - az;
    let max_var = if NEED_MAX {
        abs_max(vec![
            in bump;
            px, py, pz,
            a11, a12, a13, a21, a22, a23, a31, a32, a33,
            px_rx, py_ry, pz_rz
        ])
    } else {
        None
    };
    (Implicit3DCache { x, y, z, d }, max_var)
}

impl<'b> ImplicitPoint3D<'b> for ImplicitPointLPI<'b> {
    fn static_filter(&self) -> Option<&(Implicit3DCache<f64>, f64)> {
        if self.ss_filter.borrow().is_some() {
            let filter = self.ss_filter.borrow();
            if filter.as_ref().unwrap().1 == 0.0 {
                return None;
            } else {
                return Ref::leak(filter).as_ref();
            }
        } else {
            let (mut filter, max_var) = lpi_lambda::<'_, true, _, _>(
                self.p.data[0],
                self.p.data[1],
                self.p.data[2],
                self.q.data[0],
                self.q.data[1],
                self.q.data[2],
                self.r.data[0],
                self.r.data[1],
                self.r.data[2],
                self.s.data[0],
                self.s.data[1],
                self.s.data[2],
                self.t.data[0],
                self.t.data[1],
                self.t.data[2],
                abs_max,
                &self.bump,
            );
            let max_var = max_var.unwrap();
            let mut lambda_d_eps = max_var;
            lambda_d_eps *= lambda_d_eps;
            lambda_d_eps *= max_var;
            lambda_d_eps *= 4.884981308350689e-15;
            if filter.d > lambda_d_eps || filter.d < -lambda_d_eps {
                if filter.d < 0.0 {
                    filter.x = -filter.x;
                    filter.y = -filter.y;
                    filter.z = -filter.z;
                    filter.d = -filter.d;
                }
                self.ss_filter.replace(Some((filter, max_var)));
                return Ref::leak(self.ss_filter.borrow()).as_ref();
            } else {
                self.ss_filter.replace(Some((
                    Implicit3DCache {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        d: 0.0,
                    },
                    0.0,
                )));
                return None;
            }
        }
    }

    fn dynamic_filter(&self) -> Option<&Implicit3DCache<IntervalNumber>> {
        if self.d_filter.borrow().is_some() {
            let filter = self.d_filter.borrow();
            if filter.as_ref().unwrap().d.not_zero() {
                return Ref::leak(filter).as_ref();
            } else {
                return None;
            }
        } else {
            let (mut filter, _) = lpi_lambda::<false, IntervalNumber, _>(
                self.p.data[0].into(),
                self.p.data[1].into(),
                self.p.data[2].into(),
                self.q.data[0].into(),
                self.q.data[1].into(),
                self.q.data[2].into(),
                self.r.data[0].into(),
                self.r.data[1].into(),
                self.r.data[2].into(),
                self.s.data[0].into(),
                self.s.data[1].into(),
                self.s.data[2].into(),
                self.t.data[0].into(),
                self.t.data[1].into(),
                self.t.data[2].into(),
                dummy_abs_max,
                &self.bump,
            );
            if filter.d.negative() {
                filter.x.neg();
                filter.y.neg();
                filter.z.neg();
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

    fn exact(&self) -> Option<&Implicit3DCache<ExpansionNum<'b>>> {
        if self.exact.borrow().is_some() {
            let exact = self.exact.borrow();
            if exact.as_ref().unwrap().d.not_zero() {
                return Ref::leak(exact).as_ref();
            } else {
                return None;
            }
        } else {
            let (mut exact, _) = lpi_lambda::<'_, false, ExpansionNum<'_>, _>(
                vec![in self.bump; self.p.data[0]].into(),
                vec![in self.bump; self.p.data[1]].into(),
                vec![in self.bump; self.p.data[2]].into(),
                vec![in self.bump; self.q.data[0]].into(),
                vec![in self.bump; self.q.data[1]].into(),
                vec![in self.bump; self.q.data[2]].into(),
                vec![in self.bump; self.r.data[0]].into(),
                vec![in self.bump; self.r.data[1]].into(),
                vec![in self.bump; self.r.data[2]].into(),
                vec![in self.bump; self.s.data[0]].into(),
                vec![in self.bump; self.s.data[1]].into(),
                vec![in self.bump; self.s.data[2]].into(),
                vec![in self.bump; self.t.data[0]].into(),
                vec![in self.bump; self.t.data[1]].into(),
                vec![in self.bump; self.t.data[2]].into(),
                dummy_abs_max,
                &self.bump,
            );
            if exact.d.negative() {
                exact.x.neg();
                exact.y.neg();
                exact.z.neg();
                exact.d.neg();
            }
            normalize_lambda3d(&mut exact.x, &mut exact.y, &mut exact.z, &mut exact.d);

            self.exact.replace(Some(exact));
            if self.exact.borrow().as_ref().unwrap().d.not_zero() {
                return Ref::leak(self.exact.borrow()).as_ref();
            } else {
                return None;
            }
        }
    }
}

/// A point in 3D space representing the intersection of three planes.
#[derive(Clone)]
pub struct ImplicitPointTPI<'b> {
    /// plane 1
    pub v1: ExplicitPoint3D,
    pub v2: ExplicitPoint3D,
    pub v3: ExplicitPoint3D,
    /// plane 2
    pub w1: ExplicitPoint3D,
    pub w2: ExplicitPoint3D,
    pub w3: ExplicitPoint3D,
    /// plane 3
    pub u1: ExplicitPoint3D,
    pub u2: ExplicitPoint3D,
    pub u3: ExplicitPoint3D,

    ss_filter: RefCell<Option<(Implicit3DCache<f64>, f64)>>,
    d_filter: RefCell<Option<Implicit3DCache<IntervalNumber>>>,
    exact: RefCell<Option<Implicit3DCache<ExpansionNum<'b>>>>,

    bump: &'b Bump,
}

impl<'b> ImplicitPointTPI<'b> {
    pub fn new(
        v1: ExplicitPoint3D,
        v2: ExplicitPoint3D,
        v3: ExplicitPoint3D,
        w1: ExplicitPoint3D,
        w2: ExplicitPoint3D,
        w3: ExplicitPoint3D,
        u1: ExplicitPoint3D,
        u2: ExplicitPoint3D,
        u3: ExplicitPoint3D,
        bump: &'b Bump,
    ) -> Self {
        Self {
            v1,
            v2,
            v3,
            w1,
            w2,
            w3,
            u1,
            u2,
            u3,
            ss_filter: RefCell::new(None),
            d_filter: RefCell::new(None),
            exact: RefCell::new(None),
            bump,
        }
    }
}

fn tpi_lambda<'b, const NEED_MAX: bool, T: GenericNum + 'b, F>(
    ov1x: T,
    ov1y: T,
    ov1z: T,
    ov2x: T,
    ov2y: T,
    ov2z: T,
    ov3x: T,
    ov3y: T,
    ov3z: T,
    ow1x: T,
    ow1y: T,
    ow1z: T,
    ow2x: T,
    ow2y: T,
    ow2z: T,
    ow3x: T,
    ow3y: T,
    ow3z: T,
    ou1x: T,
    ou1y: T,
    ou1z: T,
    ou2x: T,
    ou2y: T,
    ou2z: T,
    ou3x: T,
    ou3y: T,
    ou3z: T,
    abs_max: F,
    bump: &'b Bump,
) -> (Implicit3DCache<T>, Option<T>)
where
    F: FnOnce(Vec<'b, T>) -> Option<T>,
{
    let v3x = ov3x - &ov2x;
    let v3y = ov3y - &ov2y;
    let v3z = ov3z - &ov2z;
    let v2x = ov2x - &ov1x;
    let v2y = ov2y - &ov1y;
    let v2z = ov2z - &ov1z;
    let w3x = ow3x - &ow2x;
    let w3y = ow3y - &ow2y;
    let w3z = ow3z - &ow2z;
    let w2x = ow2x - &ow1x;
    let w2y = ow2y - &ow1y;
    let w2z = ow2z - &ow1z;
    let u3x = ou3x - &ou2x;
    let u3y = ou3y - &ou2y;
    let u3z = ou3z - &ou2z;
    let u2x = ou2x - &ou1x;
    let u2y = ou2y - &ou1y;
    let u2z = ou2z - &ou1z;
    let nvx1 = &v2y * &v3z;
    let nvx2 = &v2z * &v3y;
    let nvx = nvx1 - nvx2;
    let nvy1 = &v3x * &v2z;
    let nvy2 = &v3z * &v2x;
    let nvy = nvy1 - nvy2;
    let nvz1 = &v2x * &v3y;
    let nvz2 = &v2y * &v3x;
    let nvz = nvz1 - nvz2;
    let nwx1 = &w2y * &w3z;
    let nwx2 = &w2z * &w3y;
    let nwx = nwx1 - nwx2;
    let nwy1 = &w3x * &w2z;
    let nwy2 = &w3z * &w2x;
    let nwy = nwy1 - nwy2;
    let nwz1 = &w2x * &w3y;
    let nwz2 = &w2y * &w3x;
    let nwz = nwz1 - nwz2;
    let nux1 = &u2y * &u3z;
    let nux2 = &u2z * &u3y;
    let nux = nux1 - nux2;
    let nuy1 = &u3x * &u2z;
    let nuy2 = &u3z * &u2x;
    let nuy = nuy1 - nuy2;
    let nuz1 = &u2x * &u3y;
    let nuz2 = &u2y * &u3x;
    let nuz = nuz1 - nuz2;
    let nwyuz1 = &nwy * &nuz;
    let nwyuz2 = &nwz * &nuy;
    let nwyuz = nwyuz1 - nwyuz2;
    let nwxuz1 = &nwx * &nuz;
    let nwxuz2 = &nwz * &nux;
    let nwxuz = nwxuz1 - nwxuz2;
    let nwxuy1 = &nwx * &nuy;
    let nwxuy2 = &nwy * &nux;
    let nwxuy = nwxuy1 - nwxuy2;
    let nvyuz1 = &nvy * &nuz;
    let nvyuz2 = &nvz * &nuy;
    let nvyuz = nvyuz1 - nvyuz2;
    let nvxuz1 = &nvx * &nuz;
    let nvxuz2 = &nvz * &nux;
    let nvxuz = nvxuz1 - nvxuz2;
    let nvxuy1 = &nvx * &nuy;
    let nvxuy2 = &nvy * &nux;
    let nvxuy = nvxuy1 - nvxuy2;
    let nvywz1 = &nvy * &nwz;
    let nvywz2 = &nvz * &nwy;
    let nvywz = nvywz1 - nvywz2;
    let nvxwz1 = &nvx * &nwz;
    let nvxwz2 = &nvz * &nwx;
    let nvxwz = nvxwz1 - nvxwz2;
    let nvxwy1 = &nvx * &nwy;
    let nvxwy2 = &nvy * &nwx;
    let nvxwy = nvxwy1 - nvxwy2;
    let p1a = &nvx * &ov1x;
    let p1b = &nvy * &ov1y;
    let p1c = &nvz * &ov1z;
    let p1ab = p1a + p1b;
    let p1 = p1ab + p1c;
    let p2a = nwx * &ow1x;
    let p2b = nwy * &ow1y;
    let p2c = nwz * &ow1z;
    let p2ab = p2a + p2b;
    let p2 = p2ab + p2c;
    let p3a = nux * &ou1x;
    let p3b = nuy * &ou1y;
    let p3c = nuz * &ou1z;
    let p3ab = p3a + p3b;
    let p3 = p3ab + p3c;
    let lxa = &p1 * &nwyuz;
    let lxb = &p3 * nvywz;
    let lxc = &p2 * nvyuz;
    let lxab = lxa + lxb;
    let x = lxab - lxc;
    let lya = &p2 * nvxuz;
    let lyb = &p3 * nvxwz;
    let lyc = &p1 * &nwxuz;
    let lybc = lyc + lyb;
    let y = lya - lybc;
    let lza = p3 * nvxwy;
    let lzb = p1 * &nwxuy;
    let lzc = p2 * nvxuy;
    let lzab = lza + lzb;
    let z = lzab - lzc;
    let da = nvx * nwyuz;
    let db = nvz * nwxuy;
    let dc = nvy * nwxuz;
    let dab = da + db;
    let d = dab - dc;

    let max_var = if NEED_MAX {
        abs_max(vec![
            in bump;
            ov1x, ov1y, ov1z,
            ow1x, ow1y, ow1z,
            ou1x, ou1y, ou1z,
            v3x, v3y, v3z, v2x, v2y, v2z,
            w3x, w3y, w3z, w2x, w2y, w2z,
            u3x, u3y, u3z, u2x, u2y, u2z,
        ])
    } else {
        None
    };
    (Implicit3DCache { x, y, z, d }, max_var)
}

impl<'b> ImplicitPoint3D<'b> for ImplicitPointTPI<'b> {
    fn static_filter(&self) -> Option<&(Implicit3DCache<f64>, f64)> {
        if self.ss_filter.borrow().is_some() {
            let filter = self.ss_filter.borrow();
            if filter.as_ref().unwrap().1 == 0.0 {
                return None;
            } else {
                return Ref::leak(filter).as_ref();
            }
        } else {
            let (mut filter, max_var) = tpi_lambda::<'_, true, _, _>(
                self.v1.data[0],
                self.v1.data[1],
                self.v1.data[2],
                self.v2.data[0],
                self.v2.data[1],
                self.v2.data[2],
                self.v3.data[0],
                self.v3.data[1],
                self.v3.data[2],
                self.w1.data[0],
                self.w1.data[1],
                self.w1.data[2],
                self.w2.data[0],
                self.w2.data[1],
                self.w2.data[2],
                self.w3.data[0],
                self.w3.data[1],
                self.w3.data[2],
                self.u1.data[0],
                self.u1.data[1],
                self.u1.data[2],
                self.u2.data[0],
                self.u2.data[1],
                self.u2.data[2],
                self.u3.data[0],
                self.u3.data[1],
                self.u3.data[2],
                abs_max,
                self.bump,
            );
            let max_var = max_var.unwrap();
            let mut lambda_d_eps = max_var;
            lambda_d_eps *= lambda_d_eps;
            lambda_d_eps *= lambda_d_eps;
            lambda_d_eps *= max_var;
            lambda_d_eps *= max_var;
            lambda_d_eps *= 8.704148513061234e-14;
            if filter.d > lambda_d_eps || filter.d < -lambda_d_eps {
                if filter.d < 0.0 {
                    filter.x = -filter.x;
                    filter.y = -filter.y;
                    filter.z = -filter.z;
                    filter.d = -filter.d;
                }
                self.ss_filter.replace(Some((filter, max_var)));
                return Ref::leak(self.ss_filter.borrow()).as_ref();
            } else {
                self.ss_filter.replace(Some((
                    Implicit3DCache {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        d: 0.0,
                    },
                    0.0,
                )));
                return None;
            }
        }
    }

    fn dynamic_filter(&self) -> Option<&Implicit3DCache<IntervalNumber>> {
        if self.d_filter.borrow().is_some() {
            let filter = self.d_filter.borrow();
            if filter.as_ref().unwrap().d.not_zero() {
                return Ref::leak(filter).as_ref();
            } else {
                return None;
            }
        } else {
            let (mut filter, _) = tpi_lambda::<false, IntervalNumber, _>(
                self.v1.data[0].into(),
                self.v1.data[1].into(),
                self.v1.data[2].into(),
                self.v2.data[0].into(),
                self.v2.data[1].into(),
                self.v2.data[2].into(),
                self.v3.data[0].into(),
                self.v3.data[1].into(),
                self.v3.data[2].into(),
                self.w1.data[0].into(),
                self.w1.data[1].into(),
                self.w1.data[2].into(),
                self.w2.data[0].into(),
                self.w2.data[1].into(),
                self.w2.data[2].into(),
                self.w3.data[0].into(),
                self.w3.data[1].into(),
                self.w3.data[2].into(),
                self.u1.data[0].into(),
                self.u1.data[1].into(),
                self.u1.data[2].into(),
                self.u2.data[0].into(),
                self.u2.data[1].into(),
                self.u2.data[2].into(),
                self.u3.data[0].into(),
                self.u3.data[1].into(),
                self.u3.data[2].into(),
                dummy_abs_max,
                &self.bump,
            );

            if filter.d.negative() {
                filter.x.neg();
                filter.y.neg();
                filter.z.neg();
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

    fn exact(&self) -> Option<&Implicit3DCache<ExpansionNum<'b>>> {
        if self.exact.borrow().is_some() {
            let exact = self.exact.borrow();
            if exact.as_ref().unwrap().d.not_zero() {
                return Ref::leak(exact).as_ref();
            } else {
                return None;
            }
        } else {
            let (mut exact, _) = tpi_lambda::<false, ExpansionNum<'_>, _>(
                vec![in self.bump; self.v1.data[0]].into(),
                vec![in self.bump; self.v1.data[1]].into(),
                vec![in self.bump; self.v1.data[2]].into(),
                vec![in self.bump; self.v2.data[0]].into(),
                vec![in self.bump; self.v2.data[1]].into(),
                vec![in self.bump; self.v2.data[2]].into(),
                vec![in self.bump; self.v3.data[0]].into(),
                vec![in self.bump; self.v3.data[1]].into(),
                vec![in self.bump; self.v3.data[2]].into(),
                vec![in self.bump; self.w1.data[0]].into(),
                vec![in self.bump; self.w1.data[1]].into(),
                vec![in self.bump; self.w1.data[2]].into(),
                vec![in self.bump; self.w2.data[0]].into(),
                vec![in self.bump; self.w2.data[1]].into(),
                vec![in self.bump; self.w2.data[2]].into(),
                vec![in self.bump; self.w3.data[0]].into(),
                vec![in self.bump; self.w3.data[1]].into(),
                vec![in self.bump; self.w3.data[2]].into(),
                vec![in self.bump; self.u1.data[0]].into(),
                vec![in self.bump; self.u1.data[1]].into(),
                vec![in self.bump; self.u1.data[2]].into(),
                vec![in self.bump; self.u2.data[0]].into(),
                vec![in self.bump; self.u2.data[1]].into(),
                vec![in self.bump; self.u2.data[2]].into(),
                vec![in self.bump; self.u3.data[0]].into(),
                vec![in self.bump; self.u3.data[1]].into(),
                vec![in self.bump; self.u3.data[2]].into(),
                dummy_abs_max,
                &self.bump,
            );
            if exact.d.negative() {
                exact.x.neg();
                exact.y.neg();
                exact.z.neg();
                exact.d.neg();
            }
            normalize_lambda3d(&mut exact.x, &mut exact.y, &mut exact.z, &mut exact.d);

            self.exact.replace(Some(exact));
            if self.exact.borrow().as_ref().unwrap().d.not_zero() {
                return Ref::leak(self.exact.borrow()).as_ref();
            } else {
                return None;
            }
        }
    }
}

#[derive(Clone)]
pub enum Point3D<'b> {
    Explicit(ExplicitPoint3D),
    LPI(ImplicitPointLPI<'b>),
    TPI(ImplicitPointTPI<'b>),
}

impl<'b> Point3D<'b> {
    #[inline(always)]
    pub fn explicit(&self) -> Option<&ExplicitPoint3D> {
        if let Point3D::Explicit(p) = self {
            Some(p)
        } else {
            None
        }
    }
}

#[test]
fn test_get_filter() {
    let p = ExplicitPoint3D {
        data: [0.0, 0.0, -1.0],
    };
    let q = ExplicitPoint3D {
        data: [0.0, 0.0, 1.0],
    };
    let r = ExplicitPoint3D {
        data: [-1.0, -1.0, 0.0],
    };
    let s = ExplicitPoint3D {
        data: [1.0, -1.0, 0.0],
    };
    let t = ExplicitPoint3D {
        data: [0.0, 1.0, 0.0],
    };
    let bump = Bump::new();
    let point = ImplicitPointLPI::new(p, q, r, s, t, &bump);
    assert!(point.static_filter().is_some());
}
