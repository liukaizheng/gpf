mod cylinder;
mod plane;
mod sphere;

pub use self::cylinder::Cylinder;
pub use self::plane::Plane;
pub use self::sphere::Sphere;

pub trait Surface {
    fn eval(&self, p: &[f64]) -> [f64; 4];
    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2];
}

pub enum Surf {
    Plane(Plane),
    Cylinder(Cylinder),
    Sphere(Sphere),
}

impl Surface for Surf {
    fn eval(&self, p: &[f64]) -> [f64; 4] {
        match self {
            Surf::Plane(surf) => surf.eval(p),
            Surf::Cylinder(surf) => surf.eval(p),
            Surf::Sphere(surf) => surf.eval(p),
        }
    }

    fn uv(&self, pt: &[f64], ref_pt: Option<&[f64]>) -> [f64; 2] {
        match self {
            Surf::Plane(plane) => plane.uv(pt, ref_pt),
            Surf::Cylinder(cylinder) => cylinder.uv(pt, ref_pt),
            Surf::Sphere(sphere) => sphere.uv(pt, ref_pt),
        }
    }
}

#[inline]
fn adjust_angle_to_reference(angle: f64, reference: f64) -> f64 {
    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
    let delta = angle - reference;
    angle - (delta / TWO_PI).round() * TWO_PI
}
