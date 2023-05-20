#![feature(test)]

use bumpalo::{collections::Vec, Bump};
use gpf::triangle::{triangulate, tetrahedralize};
use rand::{distributions::Uniform, rngs::SmallRng, Rng, SeedableRng};
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
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

fn read_points<'b>(name: &str, bump: &'b Bump) -> Vec<'b, f64> {
    let f = File::open(name).unwrap();

    let reader = BufReader::new(f);
    let mut points = Vec::new_in(bump);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut parts = line.split_whitespace();

        let x = parts.next().unwrap().parse::<f64>().unwrap();
        let y = parts.next().unwrap().parse::<f64>().unwrap();
        let z = parts.next().unwrap().parse::<f64>().unwrap();

        points.push(x);
        points.push(y);
        points.push(z);
    }
    points
}

#[test]
fn test_tetrahedralize() {
    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new_inclusive(-1.0, 1.0);
    let n_points = 10_0000;
    let bump = Bump::new();
    // let points: Vec<'_, f64> = read_points("123.xyz", &bump);
    let points = Vec::from_iter_in(rng.sample_iter(uniform).take(n_points * 3), &bump);
    let tets = tetrahedralize(&points, &bump);
    assert!(tets.tets.len() > 0);
}