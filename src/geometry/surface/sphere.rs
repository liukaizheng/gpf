use crate::math::{norm, normalize};

use super::{Surface, adjust_angle_to_reference};

pub struct Sphere {
    o: [f64; 3],
    r: f64,
}

impl Sphere {
    pub fn new(ox: f64, oy: f64, oz: f64, r: f64) -> Sphere {
        Sphere { o: [ox, oy, oz], r }
    }
}

impl Surface for Sphere {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        let d = [p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2]];
        let l = norm(&d);
        let val = l - self.r;
        [val, d[0] / l, d[1] / l, d[2] / l]
    }

    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2] {
        let mut d = [pt[0] - self.o[0], pt[1] - self.o[1], pt[2] - self.o[2]];
        normalize::<3>(&mut d);
        let mut uv = [d[1].atan2(d[0]), d[2].acos()];
        if let Some(ref_pt) = ref_pt {
            uv[0] = adjust_angle_to_reference(uv[0], ref_pt[0]);
        }
        uv
    }

    fn point(&self, uv: &[f64]) -> [f64; 3] {
        let cos_phi = uv[0].cos();
        let sin_phi = uv[0].sin();
        let cos_theta = uv[1].cos();
        let sin_theta = uv[1].sin();
        [
            self.o[0] + self.r * sin_phi * cos_theta,
            self.o[1] + self.r * sin_phi * sin_theta,
            self.o[2] + self.r * cos_phi,
        ]
    }
}
