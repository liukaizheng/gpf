use std::alloc::Allocator;

use super::Surface;

pub mod arc;
pub mod polyline;
pub mod segment;

fn comp_uv<'a, A: Allocator>(
    mut iter: impl Iterator<Item = &'a [f64]>,
    ref_pt: Option<&[f64]>,
    comp_uv_on_surf: impl Fn(&[f64], Option<&[f64]>) -> [f64; 2],
    alloc: A,
) -> Vec<f64, A> {
    let mut result = Vec::with_capacity_in(iter.size_hint().0, alloc);
    result.extend(comp_uv_on_surf(iter.next().unwrap(), ref_pt));
    for pt in iter {
        result.extend(comp_uv_on_surf(pt, Some(&result[result.len() - 2..])));
    }
    result
}

pub trait Curve {
    /// Reverse the curve.
    fn reversed(&self) -> Self;

    /// Discretize the curve into a sequence of 3D points.
    fn discrete<A: Allocator>(&self, alloc: A) -> Vec<f64, A>;

    /// Firstly discretize the curve into a sequence of 3D points,
    /// then compute the UV coordinates on the surface.
    fn approx_on_surf<S: Surface + ?Sized, A: Allocator + Copy>(
        &self,
        surf: &S,
        reversed: bool,
        ref_pt: Option<&[f64]>,
        alloc: A,
    ) -> Vec<f64, A> {
        let approx_points = self.discrete(alloc);

        if reversed {
            comp_uv(
                approx_points.chunks(3).rev(),
                ref_pt,
                |pt, start| surf.uv(pt, start),
                alloc,
            )
        } else {
            comp_uv(
                approx_points.chunks(3),
                ref_pt,
                |pt, start| surf.uv(pt, start),
                alloc,
            )
        }
    }
}

#[derive(Clone, Debug)]
pub enum Crv {
    Arc(arc::Arc),
    Segment(segment::Segment),
    Polyline(polyline::Polyline),
    None,
}

impl Default for Crv {
    fn default() -> Self {
        Crv::None
    }
}

impl Crv {
    #[inline]
    pub fn is_segment_or_none(&self) -> bool {
        match self {
            Crv::Segment(_) | Crv::None => true,
            _ => false,
        }
    }

    #[inline]
    pub fn is_none(&self) -> bool {
        match self {
            Crv::None => true,
            _ => false,
        }
    }
}

impl Curve for Crv {
    fn reversed(&self) -> Self {
        match self {
            Crv::Arc(arc) => Crv::Arc(arc.reversed()),
            Crv::Segment(segment) => Crv::Segment(segment.reversed()),
            Crv::Polyline(polyline) => Crv::Polyline(polyline.reversed()),
            Crv::None => Crv::None,
        }
    }

    fn discrete<A: Allocator>(&self, alloc: A) -> Vec<f64, A> {
        match self {
            Crv::Arc(arc) => arc.discrete(alloc),
            Crv::Segment(seg) => seg.discrete(alloc),
            Crv::Polyline(poly) => poly.discrete(alloc),
            Crv::None => unreachable!("None curve cannot be discretized"),
        }
    }
}
