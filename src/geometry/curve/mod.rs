use std::alloc::Allocator;

use super::{surface::Surf, Surface};

pub mod arc;
pub mod polyline;
pub mod segment;

pub trait Curve {
    fn discrete<A: Allocator>(&self, alloc: A) -> Vec<f64, A>;
    fn approx_on_surface<A: Allocator + Copy>(
        &self,
        surf: Surf,
        ref_pt: Option<&[f64]>,
        alloc: A,
    ) -> Vec<f64, A> {
        let approx_points = self.discrete(alloc);
        let n_points = approx_points.len() / 3;
        let mut uv_points = Vec::with_capacity_in(n_points << 1, alloc);
        uv_points.extend_from_slice(&surf.uv(&approx_points, ref_pt));
        for i in 1..n_points {
            let uv = surf.uv(&approx_points, Some(&uv_points[(i - 1) * 2..]));
            uv_points.extend_from_slice(&uv);
        }
        uv_points
    }
}
