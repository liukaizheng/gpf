use gpf::mesh::{Face, Halfedge, Mesh, SurfaceMesh, Vertex};

#[test]
fn build_nomanifold_mesh() {
    let bump = bumpalo::Bump::new();
    let mesh = SurfaceMesh::from(bumpalo::vec![in &bump;
        bumpalo::vec![in &bump;0, 1, 2],
        bumpalo::vec![in &bump;0, 2, 3],
        bumpalo::vec![in &bump;0, 3, 1],
        bumpalo::vec![in &bump;0, 4, 5],
        bumpalo::vec![in &bump;0, 5, 6],
        bumpalo::vec![in &bump;0, 6, 4],
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

#[test]
fn split_edge_and_face() {
    use bumpalo::collections::Vec;
    let bump = bumpalo::Bump::new();
    let mut mesh = SurfaceMesh::from(bumpalo::vec![in &bump;
        bumpalo::vec![in &bump;0, 1, 2],
        bumpalo::vec![in &bump;0, 1, 3],
        bumpalo::vec![in &bump;0, 1, 4],
        bumpalo::vec![in &bump;0, 1, 5],
    ]);
    // split edge
    {
        mesh.split_edge(0.into());
        for fid in mesh.faces() {
            let f_verts = Vec::from_iter_in(
                mesh.face(fid).halfedges().map(|hid| mesh.he_vertex(hid)),
                &bump,
            );
            assert_eq!(f_verts.len(), 4);
        }

        let new_v = 6.into();
        let mut n_new_halfedges = 0;
        for hid in mesh.vertex(new_v).incoming_halfedge() {
            assert_eq!(mesh.he_tip_vertex(hid), new_v);
            n_new_halfedges += 1;
        }
        assert_eq!(n_new_halfedges, 4);
        n_new_halfedges = 0;
        for hid in mesh.vertex(new_v).outgoing_halfedge() {
            assert_eq!(mesh.he_vertex(hid), new_v);
            n_new_halfedges += 1;
        }
        assert_eq!(n_new_halfedges, 4);
    }
    // split face
    {
        let new_hid = mesh.split_face(0.into(), 2.into(), 6.into());
        assert_eq!(
            mesh.halfedge(new_hid).face().unwrap().halfedges().count(),
            3,
        );
    }
}
