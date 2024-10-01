#![feature(allocator_api)]
use gpf::triangle::convex_2;
use itertools::Itertools;
use rand::{distributions::Uniform, rngs::SmallRng, Rng, SeedableRng};
use std::io::Write; // Add this line to import the Write trait

fn write_obj(name: &str, points: &[f64], hull: &[usize]) {
    let mut file = std::fs::File::create(name).unwrap();
    for i in 0..points.len() / 2 {
        writeln!(file, "v {} {} 0", points[i * 2], points[i * 2 + 1]).unwrap();
    }
    for (va, vb) in hull.iter().circular_tuple_windows() {
        writeln!(file, "l {} {}", va + 1, vb + 1).unwrap();
    }
}

/*#[test]
fn test_convex_2() {
    let rng = SmallRng::seed_from_u64(5489);
    let uniform = Uniform::new_inclusive(-1.0, 1.0);
    let n_points = 100_0000;
    let points = Vec::from_iter(rng.sample_iter(uniform).take(n_points * 2));

    let hull = convex_2(&points, std::alloc::Global);
    // write_obj("convex_2.obj", &points, &hull);
}

#[test]
fn test_convex_colinear() {
    let n_points = 11;
    let points: Vec<_> = (0..n_points).flat_map(|i| [i as f64, 0.0]).collect();

    let hull = convex_2(&points, std::alloc::Global);
    println!("{:?}", hull);
}*/

#[test]
fn test_convex_2_bug() {
    let points = vec![
        0.0066809675417891512,
        0.023437500000000222,
        0.0601942089234202,
        0.075000000000000178,
        -0.015322529128635853,
        0.023437500000000222,
        0.045154762758645717,
        -0.028124999999999845,
        0.02429332798771628,
        0.040625000000000203,
        -0.00094097769405939932,
        0.023437500000000222,
        0.019504265622477833,
        0.0062500000000001998,
        0.035016016201510486,
        0.057812500000000197,
        0.05492328431769869,
        0.040625000000000175,
        0.042141470936408021,
        0.057812500000000197,
        0.0044920738228399829,
        0.0062500000000001998,
        -0.0082813076490849472,
        0.023437500000000222,
        0.0098433561664043774,
        0.040625000000000203,
        0.032328907715148084,
        -0.010937499999999822,
        0.049906750196458058,
        0.0062500000000001582,
        0.024669364366963355,
        -0.010937499999999822,
        0.029708304175730393,
        0.023437500000000184,
        0.01185754801743806,
        0.0062500000000001998,
        0.037102845156667563,
        0.023437500000000184,
        0.016925863597794785,
        0.04062500000000021,
    ];

    let hull = convex_2(&points, std::alloc::Global);
    write_obj("hull.obj", &points, &hull);
}
