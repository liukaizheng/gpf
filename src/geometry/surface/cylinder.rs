use crate::math::{dot, norm};

use super::{adjust_angle_to_reference, Plane, Surface};

pub struct Cylinder {
    plane: Plane,
    r: f64,
}

impl Cylinder {
    pub fn new(ox: f64, oy: f64, oz: f64, dzx: f64, dzy: f64, dzz: f64, r: f64) -> Cylinder {
        Cylinder {
            plane: Plane::new(ox, oy, oz, dzx, dzy, dzz),
            r,
        }
    }

    #[inline]
    pub fn o(&self) -> &[f64; 3] {
        &self.plane.o
    }

    #[inline]
    pub fn dx(&self) -> &[f64; 3] {
        &self.plane.dx
    }

    #[inline]
    pub fn dy(&self) -> &[f64; 3] {
        &self.plane.dy
    }

    #[inline]
    pub fn dz(&self) -> &[f64; 3] {
        &self.plane.dz
    }
}

impl Surface for Cylinder {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        let o = self.o();
        let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
        let dz = self.dz();
        let h = dot(&d, dz);
        let d = [d[0] - h * dz[0], d[1] - h * dz[1], d[2] - h * dz[2]];
        let l = norm(&d);
        if l == 0.0 {
            [-self.r, 0.0, 0.0, 0.0]
        } else {
            [l - self.r, d[0] / l, d[1] / l, d[2] / l]
        }
    }

    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2] {
        let o = self.o();
        let vec = [pt[0] - o[0], pt[1] - o[1], pt[2] - o[2]];
        let dx = dot(self.dx(), &vec);
        let dy = dot(self.dy(), &vec);
        let mut uv = [dy.atan2(dx), dot(&vec, self.dz())];
        if let Some(ref_pt) = ref_pt {
            uv[0] = adjust_angle_to_reference(uv[0], ref_pt[0]);
        }
        uv
    }
}
