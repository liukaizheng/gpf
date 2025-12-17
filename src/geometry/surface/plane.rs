use std::alloc::Allocator;

use itertools::Itertools;

use crate::{
    geometry::{is_parallel, BBox, Crv, Curve},
    math::{add_with_coeff, cross, dot, normalize},
};

use super::{BOX_SCALE_FACTOR, Surface};

#[derive(Clone, Debug)]
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

    #[inline]
    pub fn from_x_y(o: [f64; 3], dx: [f64; 3], dy: [f64; 3]) -> Plane {
        let dz = cross(&dx, &dy);
        Plane { o, dx, dy, dz }
    }

    #[inline]
    pub fn reversed(&self) -> Plane {
        let dy = self.dy.map(|x| -x);
        Plane::from_x_y(self.o, self.dx, dy)
    }

    pub fn compute_box<'a, T: IntoIterator<Item = &'a Crv>, A: Allocator + Copy>(
        &self,
        outer_wire_edge_curves: T,
        alloc: A,
    ) -> BBox {
        let mut bbox = BBox::default();
        for crv in outer_wire_edge_curves {
            match crv {
                Crv::Segment(seg) => {
                    bbox.extend(&seg.start);
                    bbox.extend(&seg.end);
                }
                crv => {
                    bbox.merge(
                        &BBox::from_iter(crv.discrete(alloc).chunks(3))
                            .scaled(BOX_SCALE_FACTOR),
                    );
                }
            }
        }
        bbox
    }
}

impl Surface for Plane {
    fn dist(&self, p: &[f64]) -> f64 {
        let d = [p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2]];
        dot(&d, &self.dz)
    }
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

    fn equal(&self, other: &Self, tol: &crate::Tolerance) -> bool {
        if self.dist(&other.o).abs() < tol.dist_tol && other.dist(&self.o).abs() < tol.dist_tol {
            is_parallel(&self.dz, &other.dz, tol.cos_tol)
        } else {
            false
        }
    }
}
