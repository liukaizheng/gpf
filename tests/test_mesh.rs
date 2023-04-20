use gpf::mesh::{Mesh, SurfaceMesh, Vertex};

#[test]
fn build_nomanifold_mesh() {
    let mesh = SurfaceMesh::from(vec![
        vec![0, 1, 2],
        vec![0, 2, 3],
        vec![0, 3, 1],
        vec![0, 4, 5],
        vec![0, 5, 6],
        vec![0, 6, 4],
    ]);

    let base_vertices = vec![1, 2, 3, 4, 5, 6];

    // vertices of incoming halfedges
    let mut incoming_vertices = mesh
        .vertex(0.into())
        .incoming_halfedge()
        .map(|hid| *mesh.he_vertex(hid))
        .collect::<Vec<_>>();
    incoming_vertices.sort();
    assert_eq!(base_vertices, incoming_vertices);

    // vertices of outgoing halfedges
    let mut outgoing_vertices = mesh
        .vertex(0.into())
        .outgoing_halfedge()
        .map(|hid| *mesh.he_tip_vertex(hid))
        .collect::<Vec<_>>();
    outgoing_vertices.sort();
    assert_eq!(outgoing_vertices, base_vertices);

    // adjacent vertices
    let mut adjacent_vertices = mesh
        .vertex(0.into())
        .adjacent_vertices()
        .map(|vid| *vid)
        .collect::<Vec<_>>();
    adjacent_vertices.sort();
    assert_eq!(adjacent_vertices, base_vertices);
}
