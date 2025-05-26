use std::collections::BTreeMap;
use std::ops::Bound::*;

use tinyvec::{TinyVec, tiny_vec};

use crate::Tolerance;

use super::Surface;
struct Num(f64);

impl PartialEq for Num {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for Num {}

impl PartialOrd for Num {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for Num {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl From<f64> for Num {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

pub struct UniqueSurface<S: Surface> {
    surfaces: BTreeMap<Num, (S, TinyVec<[usize; 1]>)>,
}

impl<S: Surface + Clone> UniqueSurface<S> {
    pub fn new() -> Self {
        UniqueSurface {
            surfaces: BTreeMap::new(),
        }
    }

    pub fn add_surf(&mut self, surf: &S, index: usize, tol: &Tolerance) {
        const ZERO: [f64; 3] = [0.0, 0.0, 0.0];
        let dist = surf.dist(&ZERO).abs();
        for (_, (s, indices)) in self.surfaces.range_mut((
            Excluded(Num(dist - tol.dist_tol)),
            Excluded((dist + tol.dist_tol).into()),
        )) {
            if s.equal(surf, tol) {
                indices.push(index);
                return;
            }
        }

        self.surfaces
            .insert(dist.into(), (surf.clone(), tiny_vec!([usize; 1] => index)));
    }

    #[inline]
    pub fn surfaces(self) -> impl IntoIterator<Item = (S, TinyVec<[usize; 1]>)> {
        self.surfaces.into_values()
    }
}

#[cfg(test)]
mod tests {
    use crate::{Tolerance, geometry::Plane, math::normalize};

    use super::UniqueSurface;

    #[test]
    fn test_map() {
        let tol = Tolerance::new(1e-3, 1e-3);
        let plane1 = Plane::new(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let mut dz = [1e-6, 0.0, 1.0 - 1e-6];
        normalize::<3>(&mut dz);
        let plane2 = Plane::new(1.0, 0.0, 0.0, dz[0], dz[1], dz[2]);
        let mut unique_surfaces = UniqueSurface::<Plane>::new();
        unique_surfaces.add_surf(&plane1, 0, &tol);
        unique_surfaces.add_surf(&plane2, 1, &tol);
        debug_assert_eq!(unique_surfaces.surfaces.len(), 1);
    }
}
