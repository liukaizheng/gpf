use std::alloc::Allocator;

use super::Curve;

#[derive(Clone, Debug)]
pub struct Segment {
    start: [f64; 3],
    end: [f64; 3],
}

impl Segment {
    pub fn new(start: [f64; 3], end: [f64; 3]) -> Self {
        Segment { start, end }
    }
}

impl Curve for Segment {
    #[inline]
    fn reversed(&self) -> Self {
        Segment {
            start: self.end,
            end: self.start,
        }
    }

    fn discrete<A: Allocator>(&self, alloc: A) -> Vec<f64, A> {
        let mut points = Vec::new_in(alloc);
        points.extend_from_slice(&self.start);
        points.extend_from_slice(&self.end);
        points
    }
}
