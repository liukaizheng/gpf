use std::cmp::Ordering;
use std::ops::Bound::*;

use tinyvec::{TinyVec, tiny_vec};

use crate::Tolerance;
use crate::utils::RBTree;

use super::Surface;

struct DistAndSurf<S: Surface> {
    d: f64,
    surf: S,
    surf_indices: TinyVec<[usize; 1]>,
}

impl<S: Surface> PartialEq for DistAndSurf<S> {
    fn eq(&self, other: &Self) -> bool {
        self.d == other.d
    }
}

impl<S: Surface> PartialOrd for DistAndSurf<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.d.partial_cmp(&other.d)
    }
}

impl<S: Surface> PartialEq<f64> for DistAndSurf<S> {
    fn eq(&self, other: &f64) -> bool {
        self.d == *other
    }
}

impl<S: Surface> PartialOrd<f64> for DistAndSurf<S> {
    fn partial_cmp(&self, other: &f64) -> Option<Ordering> {
        self.d.partial_cmp(other)
    }
}

impl<S: Surface> PartialEq<DistAndSurf<S>> for f64 {
    fn eq(&self, other: &DistAndSurf<S>) -> bool {
        *self == other.d
    }
}

impl<S: Surface> PartialOrd<DistAndSurf<S>> for f64 {
    fn partial_cmp(&self, other: &DistAndSurf<S>) -> Option<Ordering> {
        self.partial_cmp(&other.d)
    }
}

pub struct UniqueSurface<S: Surface> {
    inner: RBTree<DistAndSurf<S>, std::alloc::Global>,
}

impl<S: Surface + Clone> UniqueSurface<S> {
    pub fn new() -> Self {
        UniqueSurface {
            inner: RBTree::new(std::alloc::Global, false),
        }
    }

    pub fn add_surf(&mut self, surf: &S, index: usize, tol: &Tolerance) {
        const ZERO: [f64; 3] = [0.0, 0.0, 0.0];
        let dist = surf.dist(&ZERO).abs();
        for entry in self
            .inner
            .range_mut((Excluded(dist - tol.dist_tol), Excluded(dist + tol.dist_tol)))
        {
            if entry.surf.equal(surf, tol) {
                entry.surf_indices.push(index);
                return;
            }
        }

        self.inner.insert(DistAndSurf {
            d: dist,
            surf: surf.clone(),
            surf_indices: tiny_vec!([usize; 1] => index),
        });
    }

    #[inline]
    pub fn into_surfaces(self) -> impl IntoIterator<Item = (S, TinyVec<[usize; 1]>)> {
        self.inner.into_iter().map(|e| (e.surf, e.surf_indices))
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
        let surfaces = unique_surfaces
            .into_surfaces()
            .into_iter()
            .collect::<Vec<_>>();
        debug_assert_eq!(surfaces.len(), 1);
    }
}
