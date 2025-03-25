use itertools::Itertools;

use crate::math::{add_with_coeff, cross, dot, normalize};

use super::Surface;

pub struct Plane {
    pub o: [f64; 3],
    pub dx: [f64; 3],
    pub dy: [f64; 3],
    pub dz: [f64; 3],
}

impl Plane {
    pub fn new(ox: f64, oy: f64, oz: f64, dzx: f64, dzy: f64, dzz: f64) -> Plane {
        let dz = [dzx, dzy, dzz];
        let min_index = dz
            .iter()
            .map(|x| x.abs())
            .position_min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap();
        let mut dx = [0.0, 0.0, 0.0];
        let i = (min_index + 1) % 3;
        let j = (i + 1) % 3;
        dx[i] = dz[j];
        dx[j] = -dz[i];
        normalize::<3>(&mut dx);
        let dy = cross(&dz, &dx);
        Plane {
            o: [ox, oy, oz],
            dx,
            dy,
            dz,
        }
    }
}

impl Surface for Plane {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        let d = [p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2]];
        let dz = &self.dz;
        [dot(&d, dz), dz[0], dz[1], dz[2]]
    }

    fn uv(&self, pt: &[f64], _ref_pt: Option<&[f64]>) -> [f64; 2] {
        let d = [pt[0] - self.o[0], pt[1] - self.o[1], pt[2] - self.o[2]];
        let dx = &self.dx;
        let dy = &self.dy;
        [dot(&d, dx), dot(&d, dy)]
    }

    fn point(&self, uv: &[f64]) -> [f64; 3] {
        add_with_coeff([(&self.dx, uv[0]), (&self.dy, uv[1]), (&self.o, 1.0)])
    }
}
