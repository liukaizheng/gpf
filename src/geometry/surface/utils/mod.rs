use super::Surface;

pub struct UniqueSurface<S: Surface> {
    surfaces: Vec<(S, Vec<usize>)>,
}

impl<S: Surface> UniqueSurface<S> {
    pub fn new() -> Self {
        UniqueSurface {
            surfaces: Vec::new(),
        }
    }

    pub fn add_surf(&mut self, surf: &S) {
        const ZERO: [f64; 3] = [0.0, 0.0, 0.0];
        let dist = surf.eval(&ZERO)[0].abs();
    }
}
