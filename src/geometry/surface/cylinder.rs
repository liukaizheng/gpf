use crate::{
    geometry::{sq_dist_to_line, is_parallel},
    math::{add_with_coeff, dot, norm},
};

use super::{Plane, Surface, adjust_angle_to_reference};

#[derive(Clone, Debug)]
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
    fn dist(&self, p: &[f64]) -> f64 {
        let o = self.o();
        let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
        let dz = self.dz();
        let h = dot(&d, dz);
        let d = [d[0] - h * dz[0], d[1] - h * dz[1], d[2] - h * dz[2]];
        let l = norm(&d);
        l - self.r
    }

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

    fn point(&self, uv: &[f64]) -> [f64; 3] {
        add_with_coeff([
            (self.dx(), uv[0].cos()),
            (self.dy(), uv[0].sin()),
            (self.dz(), uv[1]),
            (self.o(), 1.0),
        ])
    }

    fn equal(&self, other: &Self, tol: &crate::Tolerance) -> bool {
        if (self.r - other.r).abs() >= tol.dist_tol {
            return false;
        }
        is_parallel(self.dz(), other.dz(), tol.cos_tol)
            && sq_dist_to_line(self.o(), other.o(), other.dz()).abs() < tol.sq_tol
            && sq_dist_to_line(other.o(), self.o(), self.dz()).abs() < tol.sq_tol
    }
}
