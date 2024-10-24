#![feature(allocator_api)]
use gpf::triangle::convex_3;
use rand::{distributions::Uniform, rngs::SmallRng, Rng, SeedableRng};

fn write_obj(name: &str, points: &[f64], triangles: &[usize]) {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(name).unwrap();
    for i in 0..points.len() / 3 {
        writeln!(
            file,
            "v {} {} {}",
            points[i * 3],
            points[i * 3 + 1],
            points[i * 3 + 2]
        )
        .unwrap();
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
fn test_convex_3() {
    let points = vec![
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.2, 0.2, 1.0, 1.0, 1.0,
    ];

    match convex_3(&points, false, std::alloc::Global) {
        gpf::triangle::Convex3Result::Dim3(hull) => {
            write_obj("hull.obj", &points, &hull);
        }
        _ => panic!("Expected Dim3"),
    }
}


#[test]
fn test_convex_3_rng() {
    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new_inclusive(-1.0, 1.0);
    let n_points = 1000 ;
    let points = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 3));

    match convex_3(&points, false, std::alloc::Global) {
        gpf::triangle::Convex3Result::Dim3(hull) => {
            write_obj("hull.obj", &points, &hull);
        }
        _ => panic!("Expected Dim3"),
    }
}
