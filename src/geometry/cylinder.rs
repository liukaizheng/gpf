use crate::math::{dot, norm, square_norm};

use super::Surface;

pub struct Cylinder {
	o: [f64; 3],
	dz: [f64; 3],
	r: f64,
}

impl Cylinder {
	pub fn new(ox: f64, oy: f64, oz: f64, dzx: f64, dzy: f64, dzz: f64, r: f64) -> Cylinder {
		Cylinder {
			o: [ox, oy, oz],
			dz: [dzx, dzy, dzz],
			r,
		}
	}
}

impl Surface for Cylinder {
	fn eval(&self, p: &[f64]) -> [f64; 4] {
		let d = [p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2]];
		let dz = self.dz;
		let h = dot(&d, &dz);
		let d = [d[0] - h * dz[0], d[1] - h * dz[1], d[2] - h * dz[2]];
		let l = norm(&d);
		if l == 0.0 {
			[-self.r, 0.0, 0.0, 0.0]
		} else {
			[l - self.r, d[0] / l, d[1] / l, d[2] / l]
		}
	}
}