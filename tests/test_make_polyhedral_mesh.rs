use bumpalo::{collections::Vec, Bump};
use gpf::polygonlization::make_polyhedra_mesh;
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