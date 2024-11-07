use std::alloc::Allocator;

use super::{abs_max, dummy_abs_max, ExpansionNum, GenericNum, IntervalNumber};

pub fn det4<A: Allocator + Copy>(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    g: f64,
    h: f64,
    i: f64,
    j: f64,
    k: f64,
    l: f64,
    m: f64,
    n: f64,
    o: f64,
    p: f64,
    alloc: A,
) -> f64 {
    {
        // static filter
        let (det, max_var) =
            det4_impl::<true, _, _>(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, abs_max);
        let mut epsilon = max_var.unwrap();

        epsilon *= epsilon;
        epsilon *= epsilon;
        epsilon *= 1.953992523340277e-14;
        if det > epsilon || det < -epsilon {
            return det;
        }
    }
    {
        // dynamic filter
        let (det, _) = det4_impl::<false, IntervalNumber, _>(
            a.into(),
            b.into(),
            c.into(),
            d.into(),
            e.into(),
            f.into(),
            g.into(),
            h.into(),
            i.into(),
            j.into(),
            k.into(),
            l.into(),
            m.into(),
            n.into(),
            o.into(),
            p.into(),
            dummy_abs_max,
        );
        if det.not_zero() {
            return det.round();
        }
    }
    {
        // exact
        let (det, _) = det4_impl::<false, ExpansionNum<A>, _>(
            [a].to_vec_in(alloc).into(),
            [b].to_vec_in(alloc).into(),
            [c].to_vec_in(alloc).into(),
            [d].to_vec_in(alloc).into(),
            [e].to_vec_in(alloc).into(),
            [f].to_vec_in(alloc).into(),
            [g].to_vec_in(alloc).into(),
            [h].to_vec_in(alloc).into(),
            [i].to_vec_in(alloc).into(),
            [j].to_vec_in(alloc).into(),
            [k].to_vec_in(alloc).into(),
            [l].to_vec_in(alloc).into(),
            [m].to_vec_in(alloc).into(),
            [n].to_vec_in(alloc).into(),
            [o].to_vec_in(alloc).into(),
            [p].to_vec_in(alloc).into(),
            dummy_abs_max,
        );
        return *det.vec.last().unwrap();
    }
}

fn det4_impl<const NEED_MAX: bool, T: GenericNum, F>(
    a: T,
    b: T,
    c: T,
    d: T,
    e: T,
    f: T,
    g: T,
    h: T,
    i: T,
    j: T,
    k: T,
    l: T,
    m: T,
    n: T,
    o: T,
    p: T,
    abs_max: F,
) -> (T, Option<T>)
where
    F: FnOnce(&[T]) -> Option<T>,
{
    let af = &a * &f;
    let be = &b * &e;
    let kp = &k * &p;
    let lo = &l * &o;
    let ce = &c * &e;
    let ag = &a * &g;
    let jp = &j * &p;
    let ln = &l * &n;
    let ah = &a * &h;
    let de = &d * &e;
    let jo = &j * &o;
    let kn = &k * &n;
    let bg = &b * &g;
    let cf = &c * &f;
    let ip = &i * &p;
    let lm = &l * &m;
    let df = &d * &f;
    let bh = &b * &h;
    let io = &i * &o;
    let km = &k * &m;
    let ch = &c * &h;
    let dg = &d * &g;
    let in_ = &i * &n;
    let jm = &j * &m;
    let d1 = af - be;
    let d2 = kp - lo;
    let d3 = ce - ag;
    let d4 = jp - ln;
    let d5 = ah - de;
    let d6 = jo - kn;
    let d7 = bg - cf;
    let d8 = ip - lm;
    let d9 = df - bh;
    let d10 = io - km;
    let d11 = ch - dg;
    let d12 = in_ - jm;
    let t1 = d1 * d2;
    let t2 = d3 * d4;
    let t3 = d5 * d6;
    let t4 = d7 * d8;
    let t5 = d9 * d10;
    let t6 = d11 * d12;
    let r12 = t1 + t2;
    let r34 = t3 + t4;
    let r56 = t5 + t6;
    let r1234 = r12 + r34;
    let r123456 = r1234 + r56;

    let max_var = if NEED_MAX {
        abs_max(&[a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p])
    } else {
        None
    };
    (r123456, max_var)
}
