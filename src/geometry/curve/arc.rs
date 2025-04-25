use std::f64::consts::PI;

use crate::geometry::Plane;

use super::Curve;

#[derive(Clone, Debug)]
pub struct Arc {
    plane: Plane,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
}

impl Arc {
    pub fn new(plane: Plane, radius: f64, start_angle: f64, end_angle: f64) -> Self {
        Arc {
            plane,
            radius,
            start_angle,
            end_angle,
        }
    }
}

impl Curve for Arc {
    #[inline]
    fn reversed(&self) -> Self {
        Arc {
            plane: self.plane.reversed(),
            radius: self.radius,
            start_angle: -self.end_angle,
            end_angle: -self.start_angle,
        }
    }

    fn discrete<A: std::alloc::Allocator>(&self, alloc: A) -> Vec<f64, A> {
        const STEP: f64 = PI / 12.0;
        let n_steps = (((self.end_angle - self.start_angle) / STEP).round() as usize).max(1);
        let step = (self.end_angle - self.start_angle) / n_steps as f64;
        let mut points = Vec::with_capacity_in((n_steps + 1) * 3, alloc);
        for i in 0..=n_steps {
            let angle = self.start_angle + i as f64 * step;
            let x_part = self.radius * angle.cos();
            let y_part = self.radius * angle.sin();
            let x = self.plane.dx.map(|e| e * x_part);
            let y = self.plane.dy.map(|e| e * y_part);
            let o = &self.plane.o;
            points.push(o[0] + x[0] + y[0]);
            points.push(o[1] + x[1] + y[1]);
            points.push(o[2] + x[2] + y[2]);
        }
        points
    }
}
