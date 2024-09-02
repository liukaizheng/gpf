mod plane;
mod cylinder;
mod bbox;

pub use plane::Plane;
pub use cylinder::Cylinder;
pub use bbox::BBox;

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