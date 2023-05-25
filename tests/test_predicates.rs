use gpf::predicates::orient2d;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

fn read_orient2d_data(name: &str) -> (Vec<[f64; 2]>, Vec<f64>) {
    let file = File::open(name).expect("should found file");
    let reader = BufReader::new(file);
    let mut points = Vec::new();
    let mut vals = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let data = line
            .split_whitespace()
            .map(|str| str.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        points.push([data[0], data[1]]);
        vals.push(data[2]);
    }
    (points, vals)
}

#[test]
fn test_orient2d() {
    let (points, vals) = read_orient2d_data("tests/data/orient2d.txt");
    let bump = bumpalo::Bump::new();
    for i in 0..1000 {
        for j in i + 1..1000 {
            for k in j + 1..1000 {
                let o = orient2d(&points[i], &points[j], &points[k], &bump);
            }
        }
    }
}
