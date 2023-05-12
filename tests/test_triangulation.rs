#![feature(test)]

use std::time::Instant;
use bumpalo::{Bump, vec, collections::Vec};
use rand::{Rng, seq::SliceRandom};
use gpf::triangle::triangulate;
extern crate test;

#[test]
fn test_triangulate() {
    let mut rng = rand::thread_rng();
    let n_points = 500_0000;
    let bump = Bump::new();
    let mut points: Vec<'_, f64> = Vec::from_iter_in((0..n_points).map(|_|  rng.gen_range(-1.0..1.0)), &bump);
    let start = Instant::now();
    triangulate(&points, &[], &bump);
    println!("Time elapsed in {:?}", start.elapsed());
}