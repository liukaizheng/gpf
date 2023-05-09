use std::{
    cell::{Ref, RefCell},
    ops::{Add, Deref, Mul, Sub},
};

use super::{ExpansionNum, GenericNum, IntervalNumber};

pub struct ExplicitPoint3D {
    pub data: [f64; 3],
}

impl Deref for ExplicitPoint3D {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

pub struct Implicit3DCache<T> {
    x: T,
    y: T,
    z: T,
    d: T,
}

pub trait ImplicitPoint3D {
    fn static_filter(&self) -> Option<&(Implicit3DCache<f64>, f64)>;
    fn dynamic_filter(&self) -> Option<&Implicit3DCache<IntervalNumber>>;
}

/// A point in 3D space representing the intersection of a line and a plane.
pub struct ImplicitPointLPI<'b> {
    /// line
    p: ExplicitPoint3D,
    q: ExplicitPoint3D,
    /// plane
    r: ExplicitPoint3D,
    s: ExplicitPoint3D,
    t: ExplicitPoint3D,

    ssfilter: RefCell<Option<(Implicit3DCache<f64>, f64)>>,
    dfilter: RefCell<Option<Implicit3DCache<IntervalNumber>>>,
    exact: RefCell<Option<Implicit3DCache<ExpansionNum<'b>>>>,
}

impl<'b> ImplicitPointLPI<'b> {
    pub fn new(
        p: ExplicitPoint3D,
        q: ExplicitPoint3D,
        r: ExplicitPoint3D,
        s: ExplicitPoint3D,
        t: ExplicitPoint3D,
    ) -> Self {
        Self {
            p,
            q,
            r,
            s,
            t,
            ssfilter: RefCell::new(None),
            dfilter: RefCell::new(None),
            exact: RefCell::new(None),
        }
    }
}

fn lpi_ssfilter(
    p: &[f64],
    q: &[f64],
    r: &[f64],
    s: &[f64],
    t: &[f64],
) -> Option<(Implicit3DCache<f64>, f64)> {
    let a11 = p[0] - q[0];
    let a12 = p[1] - q[1];
    let a13 = p[2] - q[2];
    let a21 = s[0] - r[0];
    let a22 = s[1] - r[1];
    let a23 = s[2] - r[2];
    let a31 = t[0] - r[0];
    let a32 = t[1] - r[1];
    let a33 = t[2] - r[2];
    let tv1 = a22 * a33;
    let tv2 = a23 * a32;
    let a2233 = tv1 - tv2;
    let tv3 = a21 * a33;
    let tv4 = a23 * a31;
    let a2133 = tv3 - tv4;
    let tv5 = a21 * a32;
    let tv6 = a22 * a31;
    let a2132 = tv5 - tv6;
    let tv7 = a11 * a2233;
    let tv8 = a12 * a2133;
    let tv9 = a13 * a2132;
    let tt1 = tv7 - tv8;

    let d = tt1 + tv9;

    let px_rx = p[0] - r[0];
    let py_ry = p[1] - r[1];
    let pz_rz = p[2] - r[2];
    let tt2 = py_ry * a2133;
    let tt3 = px_rx * a2233;
    let tt4 = pz_rz * a2132;
    let tt5 = tt3 + tt4;
    let n = tt5 - tt2;
    let ax = a11 * n;
    let ay = a12 * n;
    let az = a13 * n;
    let dpx = d * p[0];
    let dpy = d * p[1];
    let dpz = d * p[2];
    let x = dpx - ax;
    let y = dpy - ay;
    let z = dpz - az;

    let max_var = p[0]
        .abs()
        .max(p[1].abs())
        .max(p[2].abs())
        .max(a11.abs())
        .max(a12.abs())
        .max(a13.abs())
        .max(a21.abs())
        .max(a22.abs())
        .max(a23.abs())
        .max(a31.abs())
        .max(a32.abs())
        .max(a33.abs())
        .max(px_rx.abs())
        .max(py_ry.abs())
        .max(pz_rz.abs());
    let mut lambda_d_eps = max_var;
    lambda_d_eps *= lambda_d_eps;
    lambda_d_eps *= max_var;
    lambda_d_eps *= 4.884981308350689e-15;
    if d > lambda_d_eps || d < -lambda_d_eps {
        Some((Implicit3DCache { x, y, z, d }, max_var))
    } else {
        None
    }
}

/*fn lpi_lambda<T>(
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
) -> (Implicit3DCache<T>, Option<T>)
 where
  T: GenericNum,
  for <'a> &'a T: Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
//   for <'a> &'a T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
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
    let tt2 = py_ry * &a2133;
    let tt3 = px_rx * &a2233;
    let tt4 = pz_rz * &a2132;
    let tt5 = tt3 + tt4;
    let n = tt5 - tt2;
    let ax = &a11 * &n;
    let ay = &a12 * &n;
    let az = &a13 * &n;
    let dpx = &d * px;
    let dpy = &d * py;
    let dpz = &d * pz;
    let x = dpx - ax;
    let y = dpy - ay;
    let z = dpz - az;
    (Implicit3DCache { x, y, z, d }, None)
}*/

impl<'b> ImplicitPoint3D for ImplicitPointLPI<'b> {
    fn static_filter(&self) -> Option<&(Implicit3DCache<f64>, f64)> {
        let filter_option = Ref::leak(self.ssfilter.borrow()).as_ref();
        if let Some(filter) = filter_option {
            if filter.1 == 0.0 {
                return None;
            } else {
                return filter_option;
            }
        } else {
            let filter = lpi_ssfilter(
                &self.p.data,
                &self.q.data,
                &self.r.data,
                &self.s.data,
                &self.t.data,
            );
            if filter.is_none() {
                self.ssfilter.replace(Some((
                    Implicit3DCache {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        d: 0.0,
                    },
                    0.0,
                )));
                return None;
            } else {
                self.ssfilter.replace(filter);
                return Ref::leak(self.ssfilter.borrow()).as_ref();
            }
        }
    }

    fn dynamic_filter(&self) -> Option<&Implicit3DCache<IntervalNumber>> {
        None
    }
}

/// A point in 3D space representing the intersection of three planes.
pub struct ImplicitPointTPI<'b> {
    /// plane 1
    v1: ExplicitPoint3D,
    v2: ExplicitPoint3D,
    v3: ExplicitPoint3D,
    /// plane 2
    w1: ExplicitPoint3D,
    w2: ExplicitPoint3D,
    w3: ExplicitPoint3D,
    /// plane 3
    u1: ExplicitPoint3D,
    u2: ExplicitPoint3D,
    u3: ExplicitPoint3D,

    ssfilter: RefCell<Option<Implicit3DCache<f64>>>,
    dfilter: RefCell<Option<Implicit3DCache<IntervalNumber>>>,
    exact: RefCell<Option<Implicit3DCache<ExpansionNum<'b>>>>,
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
            ssfilter: RefCell::new(None),
            dfilter: RefCell::new(None),
            exact: RefCell::new(None),
        }
    }
}

pub enum Point3D<'b> {
    Explicit(ExplicitPoint3D),
    LPI(ImplicitPointLPI<'b>),
    TPI(ImplicitPointTPI<'b>),
}
