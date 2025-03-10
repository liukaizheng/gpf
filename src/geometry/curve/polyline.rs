use super::Curve;

pub struct Polyline {
    points: Vec<f64>,
}

impl Curve for Polyline {
    fn discrete<A: std::alloc::Allocator>(&self, alloc: A) -> Vec<f64, A> {
        let mut points = Vec::with_capacity_in(self.points.len(), alloc);
        points.extend_from_slice(&self.points);
        points
    }
}
