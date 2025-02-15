use super::Surface;

pub struct Plane {
    o: [f64; 3],
    dz: [f64; 3],
}

impl Plane {
    pub fn new(ox: f64, oy: f64, oz: f64, dzx: f64, dzy: f64, dzz: f64) -> Plane {
        Plane {
            o: [ox, oy, oz],
            dz: [dzx, dzy, dzz],
        }
    }
}

impl Surface for Plane {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        let d = [p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2]];
        let dz = &self.dz;
        let d = d[0] * dz[0] + d[1] * dz[1] + d[2] * dz[2];
        [d, dz[0], dz[1], dz[2]]
    }
}
