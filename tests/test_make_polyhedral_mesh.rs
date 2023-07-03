use gpf::{polygonlization::make_polyhedra_mesh, predicates::ExpansionNum};
use serde::Deserialize;

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct PolygonsData {
    pointsData: Vec<f64>,
    faceInShellData: Vec<usize>,
    separators: Vec<usize>,
    edgesData: Vec<usize>,
    axesData: Vec<f64>,
}

fn make_two_dim_arr<T: Clone>(arr: &[T], separators: &[usize]) -> Vec<Vec<T>> {
    Vec::from_iter(
        separators
            .windows(2)
            .into_iter()
            .map(|a| Vec::from_iter(arr[a[0]..a[1]].iter().map(|t| t.clone()))),
    )
}

fn my_sub(a: &[f64], b: &[f64]) -> ExpansionNum {
    let c = vec![0.0; 100000000];
    ExpansionNum { vec: vec![a[0] - b[0]] }
}

#[test]
fn two_add() {
    let a = ExpansionNum {
        vec: vec![0.0],
    };
    let b = ExpansionNum {
        vec: vec![0.0],
    };
    let c = my_sub(&a, &b);
    println!("{:?}", c.vec);
}

#[test]
fn two_models() {
    let f = std::fs::File::open("tests/data/solidFix.json").expect("read file");
    let reader = std::io::BufReader::new(f);
    let data = serde_json::from_reader::<_, PolygonsData>(reader).expect("failed to read");
    let edges = make_two_dim_arr(&data.edgesData, &data.separators);
    make_polyhedra_mesh(
        &data.pointsData,
        &data.axesData,
        &data.faceInShellData,
        &edges,
        1e-6,
    );
}

// #[test]
fn test_make_polyhedra_mesh() {
    let points = vec![
        0.0, 0.0, 0.0, // 0
        1.0, 0.0, 0.0, // 1
        1.0, 1.0, 0.0, // 2
        0.0, 1.0, 0.0, // 3
        0.0, 0.0, 1.0, // 4
        1.0, 0.0, 1.0, // 5
        1.0, 1.0, 1.0, // 6
        0.0, 1.0, 1.0, // 7
    ];

    let edges = vec![
        vec![0, 3, 3, 2, 2, 1, 1, 0], // bottom
        vec![4, 5, 5, 6, 6, 7, 7, 4], // top
        vec![0, 4, 4, 7, 7, 3, 3, 0], // left
        vec![1, 2, 2, 6, 6, 5, 5, 1], // right
        vec![0, 1, 1, 5, 5, 4, 4, 0], // front
    ];

    let axis = vec![
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, // bottom
        0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // top
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, // left
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, // right
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, // front
    ];

    let poly_in_shell = vec![0, 0, 0, 0, 0];
    make_polyhedra_mesh(&points, &axis, &poly_in_shell, &edges, 1e-6);
}
