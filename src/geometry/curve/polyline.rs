use itertools::Itertools;

use super::Curve;

#[derive(Clone, Debug)]
pub struct Polyline {
    points: Vec<f64>,
}

impl Polyline {
    pub fn new(points: Vec<f64>) -> Self {
        Polyline { points }
    }
}

impl Curve for Polyline {
    #[inline]
    fn reversed(&self) -> Self {
        Polyline {
            points: self
                .points
                .chunks(3)
                .rev()
                .flat_map(|p| [p[0], p[1], p[2]])
                .collect_vec(),
        }
    }

    fn discrete<A: std::alloc::Allocator>(&self, alloc: A) -> Vec<f64, A> {
        let mut points = Vec::with_capacity_in(self.points.len(), alloc);
        points.extend_from_slice(&self.points);
        points
    }
}
