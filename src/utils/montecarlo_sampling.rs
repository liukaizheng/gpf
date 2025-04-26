use std::alloc::Allocator;

use rand::{Rng, SeedableRng, distributions::Uniform, rngs::SmallRng};

use crate::{
    math::{cross, norm, sub_short},
    point,
};

pub fn montecarlo_sampling<const N: usize, A: Allocator + Copy>(
    points: &[f64],
    triangles: &[usize],
    n_sampling_points: usize,
    alloc: A,
) -> Vec<f64, A> {
    let mut result = Vec::with_capacity_in(n_sampling_points * N, alloc);
    let mut triangle_areas = Vec::with_capacity_in(triangles.len() / 3 + 1, alloc);
    triangle_areas.push(0.0);
    triangle_areas.extend(triangles.chunks(3).map(|tri| {
        let p0 = point::<N>(points, tri[0]);
        let p1 = point::<N>(points, tri[1]);
        let p2 = point::<N>(points, tri[2]);
        let v1 = sub_short::<N, _>(p1, p0);
        let v2 = sub_short::<N, _>(p2, p0);
        if N == 3 {
            norm(&cross(&v1, &v2))
        } else {
            (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        }
    }));
    {
        let sum: f64 = triangle_areas.iter().sum();
        let mut prev = 0.0;
        for area in &mut triangle_areas {
            *area = prev + *area / sum;
            prev = *area;
        }
    }

    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new(0.0, 1.0);
    result.extend(
        rng.sample_iter(uniform)
            .array_chunks::<4>()
            .take(n_sampling_points)
            .map(|nums| {
                let idx =
                    match triangle_areas.binary_search_by(|x| x.partial_cmp(&nums[0]).unwrap()) {
                        Ok(pos) => pos,
                        Err(pos) => pos - 1,
                    } * 3;
                let s = nums[1] + nums[2] + nums[3];
                let c1 = nums[1] / s;
                let c2 = nums[2] / s;
                let c3 = nums[3] / s;
                let p1 = point::<N>(points, triangles[idx]);
                let p2 = point::<N>(points, triangles[idx + 1]);
                let p3 = point::<N>(points, triangles[idx + 2]);
                let mut ret: [f64; N] = [0.0; N];
                for i in 0..N {
                    ret[i] = p1[i] * c1 + p2[i] * c2 + p3[i] * c3;
                }
                ret
            })
            .flatten(),
    );
    result
}
