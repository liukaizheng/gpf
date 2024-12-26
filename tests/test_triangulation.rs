#![feature(test)]

use bumpalo::Bump;
use gpf::triangle::{tetrahedralize, triangulate, triangulate1};
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

#[test]
fn test_simple() {
    // let points = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0];

    let bump = Bump::new();
    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new_inclusive(-1.0, 1.0);
    let n_points = 20_0000;
    let points = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 2));

    {
        let start = Instant::now();
        // let triangles = triangulate(&points, &[], &bump);
        println!("old Time elapsed in {:?}", start.elapsed());
        // write_obj(&points, &triangles, "123.obj");
    }

    {
        let start = Instant::now();
        let triangles = triangulate1(&points, &[], true, &bump);
        println!("new Time elapsed in {:?}", start.elapsed());
        // write_obj(&points, &triangles, "124.obj");
    }
}

#[test]
fn test_bug() {
    let points = [
        -0.37532790381475534,
        -0.5701866788379338,
        -0.3746855466249469,
        -0.5701664269221419,
        -0.37421247631072985,
        -0.5692436605107906,
        -0.37385918384824846,
        -0.5691439523641123,
        -0.3710912705614987,
        -0.5629674306723101,
        -0.37098555925350085,
        -0.5609913087007848,
        -0.37096292038201095,
        -0.5597273368099648,
    ];

    let bump = Bump::new();
    let triangles = triangulate1(&points, &[], false, &bump);
    write_obj(&points, &triangles, "bug.obj");
}
