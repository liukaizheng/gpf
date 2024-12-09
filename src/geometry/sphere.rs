use crate::math::norm;

use super::Surface;

pub struct Sphere {
    o: [f64; 3],
    r: f64,
}

impl Sphere {
    pub fn new(ox: f64, oy: f64, oz: f64, r: f64) -> Sphere {
        Sphere {
            o: [ox, oy, oz],
            r,
        }
    }
}

impl Surface for Sphere {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        let d = [
            p[0] - self.o[0],
            p[1] - self.o[1],
            p[2] - self.o[2],
        ];
        let l = norm(&d);
        let val = l - self.r;
        [val, d[0] / l, d[1] / l, d[2] / l]
    }
}
