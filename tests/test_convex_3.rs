#![feature(allocator_api)]
use gpf::triangle::convex_3;
use rand::{Rng, SeedableRng, distr::Uniform, rngs::SmallRng};

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
    let uniform = Uniform::new_inclusive(-1.0, 1.0).unwrap();
    let n_points = 1000;
    let points = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 3));

    match convex_3(&points, false, std::alloc::Global) {
        gpf::triangle::Convex3Result::Dim3(hull) => {
            write_obj("hull.obj", &points, &hull);
        }
        _ => panic!("Expected Dim3"),
    }
}

#[test]
fn test_convex_3_bug_1() {
    let points = [
        -0.0107421875,
        -0.001953125,
        -0.001953125,
        0.0166015625,
        0.0234375,
        0.0234375,
        0.0009765625,
        -0.001953125,
        -0.001953125,
        -0.029296875,
        -0.0283203125,
        0.0234375,
        -0.0017555245234293342,
        0.006640624999999992,
        0.0066406250000000284,
        -0.0069528950092088608,
        -0.001953125,
        -0.001953125,
        -0.016925025966998943,
        -0.010546874999999992,
        0.0066406250000000284,
        0.01116646371973545,
        0.014843750000000008,
        0.014843749999999972,
        0.0013304975486321032,
        0.006250000000000016,
        0.0234375,
        0.0075039989672254347,
        0.014843750000000008,
        0.014843749999999972,
        -0.0091472052257138645,
        -0.010546874999999992,
        0.0066406250000000284,
        -0.0029559938680561301,
        -0.001953125,
        -0.001953125,
        0.0059008667989032563,
        0.006640624999999992,
        0.0066406250000000284,
        -0.023114698291792465,
        -0.01972656250000001,
        0.014843749999999972,
        -0.014119280578955458,
        -0.011132812500000016,
        0.0234375,
        -0.019335265169942412,
        -0.01972656250000001,
        0.014843749999999972,
        -0.004097855726835231,
        -0.0022786458333333343,
        0.014973958333333334,
        -0.0130973542162615,
        -0.0107421875,
        0.006510416666666667,
        -0.0078637582113296657,
        -0.0022786458333333343,
        0.014973958333333334,
        0.0020874061046257872,
        0.006510416666666667,
        0.006510416666666667,
        0.0,
        0.0,
        0.0,
    ];
    match convex_3(&points, false, std::alloc::Global) {
        gpf::triangle::Convex3Result::Dim3(hull) => {
            write_obj("hull.obj", &points, &hull);
        }
        _ => panic!("Expected Dim3"),
    }
}
