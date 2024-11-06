mod bbox;
mod cylinder;
mod plane;

pub use bbox::BBox;
pub use cylinder::Cylinder;
pub use plane::Plane;

pub trait Surface {
    fn eval(&self, p: &[f64]) -> [f64; 4];
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
