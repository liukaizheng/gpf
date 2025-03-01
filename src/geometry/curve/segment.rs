use std::alloc::Allocator;

use super::Curve;

pub struct Segment {
    start: [f64; 3],
    end: [f64; 3],
}

impl Curve for Segment {
    fn discrete<A: Allocator>(&self, alloc: A) -> Vec<f64, A> {
        let mut points = Vec::new_in(alloc);
        points.extend_from_slice(&self.start);
        points.extend_from_slice(&self.end);
        points
    }
}
