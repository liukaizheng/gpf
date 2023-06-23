#![feature(test)]

use bumpalo::Bump;
use gpf::triangle::{tetrahedralize, triangulate};
use rand::{distributions::Uniform, rngs::SmallRng, Rng, SeedableRng};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
extern crate test;

#[allow(dead_code)]
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
    let n_points = 1_0000;
    let points: Vec<f64> = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 2));
    let indices = Vec::from_iter(0..n_points);
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
    let bump = Bump::new();
    let triangles = triangulate(&points, &[start_idx, end_idx], &bump);
    println!("Time elapsed in {:?}", start.elapsed());
    assert_eq!(triangles.len(), 54684);
    // write_obj(&points, &triangles, "test.obj");
}

#[allow(dead_code)]
fn read_points(name: &str) -> Vec<f64> {
    let f = File::open(name).unwrap();

    let reader = BufReader::new(f);
    let mut points = Vec::new();
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
    let n_points = 1_0000;
    // let points: Vec<'_, f64> = read_points("123.xyz", &bump);
    let points = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 3));
    let tets = tetrahedralize(&points);
    assert!(tets.tets.len() > 0);
}
