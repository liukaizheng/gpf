use bumpalo::Bump;
use gpf::{polygonlization::make_polyhedra_mesh, triangle::tetrahedralize};
use serde::{Deserialize, Serialize};
type BVec<'b, T> = bumpalo::collections::Vec<'b, T>;

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct PolygonsData {
    pointsData: Vec<f64>,
    faceInShellData: Vec<usize>,
    separators: Vec<usize>,
    edgesData: Vec<usize>,
    axesData: Vec<f64>,
}

fn make_two_dim_arr<'b, T: 'b + Clone>(
    arr: &[T],
    separators: &[usize],
    bump: &'b Bump,
) -> BVec<'b, BVec<'b, T>> {
    BVec::from_iter_in(
        separators
            .windows(2)
            .into_iter()
            .map(|a| BVec::from_iter_in(arr[a[0]..a[1]].iter().map(|t| t.clone()), bump)),
        bump,
    )
}

#[test]
fn two_models() {
    let f = std::fs::File::open("tests/data/solidFix.json").expect("read file");
    let reader = std::io::BufReader::new(f);
    let data = serde_json::from_reader::<_, PolygonsData>(reader).expect("failed to read");
    let bump = Bump::new();
    let edges = make_two_dim_arr(&data.edgesData, &data.separators, &bump);
    make_polyhedra_mesh(
        &data.pointsData,
        &data.axesData,
        &data.faceInShellData,
        &edges,
        &bump,
        1e-6,
    );
    let a = 2;
}

#[test]
fn test_make_polyhedra_mesh() {
    let bump = Bump::new();
    let points = bumpalo::vec![ in &bump;
        0.0, 0.0, 0.0, // 0
        1.0, 0.0, 0.0, // 1
        1.0, 1.0, 0.0, // 2
        0.0, 1.0, 0.0, // 3
        0.0, 0.0, 1.0, // 4
        1.0, 0.0, 1.0, // 5
        1.0, 1.0, 1.0, // 6
        0.0, 1.0, 1.0, // 7
    ];

    let edges = bumpalo::vec![in &bump;
        bumpalo::vec![in &bump; 0, 3, 3, 2, 2, 1, 1, 0], // bottom
        bumpalo::vec![in &bump; 4, 5, 5, 6, 6, 7, 7, 4], // top
        bumpalo::vec![in &bump; 0, 4, 4, 7, 7, 3, 3, 0], // left
        bumpalo::vec![in &bump; 1, 2, 2, 6, 6, 5, 5, 1], // right
        bumpalo::vec![in &bump; 0, 1, 1, 5, 5, 4, 4, 0], // front
    ];

    let axis = bumpalo::vec![in &bump;
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, // bottom
        0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // top
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, // left
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, // right
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, // front
    ];

    let poly_in_shell = bumpalo::vec![in &bump;
        0, 0, 0, 0,0
    ];
    make_polyhedra_mesh(&points, &axis, &poly_in_shell, &edges, &bump, 1e-6);
    let a = 2;
}
