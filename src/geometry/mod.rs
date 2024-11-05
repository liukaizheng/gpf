mod bbox;
mod cylinder;
mod plane;

pub use bbox::BBox;
pub use cylinder::Cylinder;
pub use plane::Plane;

pub trait Surface {
    fn eval(&self, p: &[f64]) -> [f64; 4];

    #[inline]
    fn eval_round(&self, p: &[f64], eps: f64) -> [f64; 4] {
        let mut res = self.eval(p);
        res[0] = (res[0] / eps).round() * eps;
        res
    }
}

pub enum Surf {
    Plane(Plane),
    Cylinder(Cylinder),
}

impl Surface for Surf {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        match self {
            Surf::Plane(surf) => surf.eval(p),
            Surf::Cylinder(surf) => surf.eval(p),
        }
    }
}
