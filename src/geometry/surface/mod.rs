mod cylinder;
mod plane;
mod sphere;

use std::alloc::Allocator;

use itertools::Itertools;

use crate::{triangle::triangulate_with_new_points, utils::montecarlo_sampling};

pub use self::cylinder::Cylinder;
pub use self::plane::Plane;
pub use self::sphere::Sphere;

use super::{BBox, Crv, Curve};

pub const BOX_SCALE_FACTOR: f64 = 1.05;

pub trait Surface {
    fn eval(&self, p: &[f64]) -> [f64; 4];
    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2];
    fn point(&self, uv: &[f64]) -> [f64; 3];

    fn compute_face_uv_loops<
        'a,
        T: IntoIterator<Item = impl IntoIterator<Item = (&'a Crv, bool)>>,
        A: Allocator + Copy,
    >(
        &self,
        face_loops: T,
        alloc: A,
    ) -> Vec<Vec<f64, A>, A> {
        let mut result = Vec::new_in(alloc);
        for face_loop in face_loops {
            let mut ref_pt = None;
            let mut loop_uv_points = Vec::new_in(alloc);
            for (crv, reversed) in face_loop {
                loop_uv_points
                    .extend_from_slice(&crv.approx_on_surf(self, reversed, ref_pt, alloc)[2..]);
                ref_pt = Some(&loop_uv_points[(loop_uv_points.len() - 2)..]);
            }
            result.push(loop_uv_points);
        }
        result
    }

    fn compute_uv_face<
        'a,
        T: IntoIterator<Item = impl IntoIterator<Item = (&'a Crv, bool)>>,
        A: Allocator + Copy,
    >(
        &self,
        face_loops: T,
        alloc: A,
    ) -> (Vec<f64, A>, Vec<usize, A>) {
        let uv_loops = self.compute_face_uv_loops(face_loops, alloc);
        let mut uv_points: Vec<f64, A> = Vec::new_in(alloc);
        for loop_uv_points in &uv_loops {
            uv_points.extend_from_slice(&loop_uv_points);
        }
        let mut segments = Vec::new_in(alloc);
        let mut start = 0;
        for loop_uv_points in uv_loops {
            let end = start + (loop_uv_points.len() >> 1);
            for (i, j) in (start..end).circular_tuple_windows() {
                segments.push(i);
                segments.push(j);
            }
            start = end;
        }
        let (new_points, triangles) =
            triangulate_with_new_points(&uv_points, &segments, true, alloc);
        uv_points.extend(new_points);
        (uv_points, triangles)
    }

    fn compute_box_from_uv_face<A: Allocator + Copy>(
        &self,
        uv_points: &[f64],
        uv_triangles: &[usize],
        alloc: A,
    ) -> BBox {
        const N_SAMPLING: usize = 200;
        let new_points = montecarlo_sampling::<2, _>(uv_points, uv_triangles, N_SAMPLING, alloc);
        BBox::from_iter(
            uv_points
                .iter()
                .map(|x| *x)
                .chain(new_points)
                .array_chunks::<2>()
                .map(|uv| self.point(&uv)),
        )
        .scaled(BOX_SCALE_FACTOR)
    }
}

#[derive(Clone, Debug)]
pub enum Surf {
    Plane(Plane),
    Cylinder(Cylinder),
    Sphere(Sphere),
}

impl Surface for Surf {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        match self {
            Surf::Plane(surf) => surf.eval(p),
            Surf::Cylinder(surf) => surf.eval(p),
            Surf::Sphere(surf) => surf.eval(p),
        }
    }

    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2] {
        match self {
            Surf::Plane(plane) => plane.uv(pt, ref_pt),
            Surf::Cylinder(cylinder) => cylinder.uv(pt, ref_pt),
            Surf::Sphere(sphere) => sphere.uv(pt, ref_pt),
        }
    }

    fn point(&self, uv: &[f64]) -> [f64; 3] {
        match self {
            Surf::Plane(plane) => plane.point(uv),
            Surf::Cylinder(cylinder) => cylinder.point(uv),
            Surf::Sphere(sphere) => sphere.point(uv),
        }
    }
}

#[inline]
fn adjust_angle_to_reference(angle: f64, reference: f64) -> f64 {
    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
    let delta = angle - reference;
    angle - (delta / TWO_PI).round() * TWO_PI
}
