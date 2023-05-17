#![feature(test)]

use bumpalo::{collections::Vec, vec, Bump};
use gpf::triangle::triangulate;
use rand::{distributions::Uniform, rngs::SmallRng, Rng, SeedableRng};
use std::io::Write;
use std::time::Instant;
extern crate test;

fn write_obj(points: &[f64], triangles: &[usize], name: &str) {
    let mut file = std::fs::File::create(name).unwrap();
    for i in 0..points.len() / 2 {
        writeln!(file, "v {} {} 0", points[i * 2], points[i * 2 + 1]).unwrap();
    }
    for i in 0..triangles.len() / 3 {
        writeln!(
            file,
            "f {} {} {}",
            triangles[i * 3] + 1,
            triangles[i * 3 + 1] + 1,
            triangles[i * 3 + 2] + 1
        )
        .unwrap();
    }
}

#[test]
fn test_triangulate() {
    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new_inclusive(-1.0, 1.0);
    let n_points = 1_00_0000;
    let bump = Bump::new();
    let points: Vec<'_, f64> =
        Vec::from_iter_in(rng.sample_iter(uniform).take(n_points * 2), &bump);
    let indices = Vec::from_iter_in(0..n_points, &bump);
    let start_idx = *indices
        .iter()
        .min_by(|&&a, &&b| {
            let i = a << 1;
            let j = b << 1;
            (points[i], points[i + 1])
                .partial_cmp(&(points[j], points[j + 1]))
                .unwrap()
        })
        .unwrap();
    let end_idx = *indices
        .iter()
        .max_by(|&&a, &&b| {
            let i = a << 1;
            let j = b << 1;
            (points[i], points[i + 1])
                .partial_cmp(&(points[j], points[j + 1]))
                .unwrap()
        })
        .unwrap();
    let start = Instant::now();
    let triangles = triangulate(&points, &[start_idx, end_idx], &bump);
    println!("Time elapsed in {:?}", start.elapsed());
    write_obj(&points, &triangles, "test.obj");
}
