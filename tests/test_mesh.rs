use gpf::mesh::{
    Edge, EdgeId, Face, FaceId, Halfedge, HalfedgeId, Mesh, SurfaceMesh, VertexId,
};

fn validate_mesh_connectivity(mesh: &SurfaceMesh) -> Result<(), String> {
    let validate_vertex = |vid: VertexId, msg: &str| {
        if vid.0 > mesh.n_vertices_capacity() || !mesh.vertex_is_valid(vid) {
            Err(format!("{} bad vertex reference: {}", msg, vid.0))
        } else {
            Ok(())
        }
    };
    let validate_halfedge = |hid: HalfedgeId, msg: &str| {
        if hid.0 > mesh.n_halfedges_capacity() || !mesh.halfedge_is_valid(hid) {
            Err(format!("{} bad halfedge reference: {}", msg, hid.0))
        } else {
            Ok(())
        }
    };
    let validate_edge = |eid: EdgeId, msg: &str| {
        if eid.0 > mesh.n_edges_capacity() || !mesh.edge_is_valid(eid) {
            Err(format!("{} bad edge reference: {}", msg, eid.0))
        } else {
            Ok(())
        }
    };
    let validate_face = |fid: FaceId, msg: &str| {
        if fid.0 > mesh.n_faces_capacity() || !mesh.face_is_valid(fid) {
            Err(format!("{} bad face reference: {}", msg, fid.0))
        } else {
            Ok(())
        }
    };

    for vid in mesh.vertices() {
        validate_halfedge(mesh.v_halfedge(*vid), "v_he: ")?;
    }

    for hid in mesh.halfedges() {
        validate_vertex(mesh.he_from(*hid), "he_vertex: ")?;

        validate_halfedge(mesh.he_next(*hid), "he_next: ")?;
        validate_halfedge(mesh.he_twin(*hid), "he_twin: ")?;
        validate_halfedge(mesh.he_next_incoming_neighbor(*hid), "next_incoming: ")?;

        validate_edge(mesh.he_edge(*hid), "he_edge: ")?;
        let fid = mesh.he_face(*hid);
        validate_face(fid, "he_face: ")?;
    }

    for eid in mesh.edges() {
        validate_halfedge(mesh.e_halfedge(*eid), "e_he: ")?;
    }

    for fid in mesh.faces() {
        validate_halfedge(mesh.f_halfedge(*fid), "e_face: ")?;
    }

    for hid in mesh.halfedges() {
        let first_he = mesh.halfedge(*hid);
        let mut curr_he = mesh.halfedge(*hid);
        curr_he = curr_he.sibling();
        let mut count = 1;
        while *curr_he != *hid {
            if *curr_he.edge() != *first_he.edge() {
                return Err(format!(
                    "halfedge sibling doesn't have edge == he.edge for he {}",
                    hid.0
                ));
            }
            if count > mesh.n_halfedges() {
                return Err(format!(
                    "halfedge sibling doesn't cycle back for he {}",
                    hid.0
                ));
            }
            curr_he = curr_he.sibling();
            count += 1;
        }
    }

    for eid in mesh.edges() {
        for hid in mesh.edge(*eid).halfedges() {
            if *eid != mesh.he_edge(*hid) {
                return Err(format!(
                    "edge.halfedge doesn't match he.edge for edge {}",
                    eid.0
                ));
            }
        }
    }

    for face in mesh.faces() {
        let he = face.halfedge();
        let fid = *face;
        if mesh.he_face(*he) != fid {
            return Err(format!("f.he().face() is not f for face {}", fid.0));
        }

        let mut curr_he = he.next();
        let mut count = 1;
        while *curr_he != *he {
            if mesh.he_face(*curr_he) != fid {
                return Err(format!("face.he doesn't match he.face for face {}", fid.0));
            }
            if count > mesh.n_halfedges() {
                return Err(format!(
                    "halfedge next doesn't cycle back for face {}",
                    fid.0
                ));
            }
            curr_he = curr_he.next();
            count += 1;
        }
        if count < 2 {
            return Err(format!("face {} degree < 2", fid.0));
        }
    }

    for hid in mesh.halfedges() {
        let hid = *hid;
        let tip = mesh.he_to(hid);
        if *mesh
            .halfedge(mesh.he_next_incoming_neighbor(hid))
            .to()
            != tip
        {
            return Err(format!("next incoming he is not to same vert"));
        }
    }

    let bump = bumpalo::Bump::new();
    let mut v_in_count = bumpalo::vec![in &bump; 0; mesh.n_vertices_capacity()];
    let mut v_out_count = bumpalo::vec![in &bump; 0; mesh.n_vertices_capacity()];
    for hid in mesh.halfedges() {
        v_out_count[mesh.he_from(*hid)] += 1;
        v_in_count[mesh.he_to(*hid)] += 1;
    }
    for vid in mesh.vertices() {
        let vid = *vid;
        let mut count = 0;
        for _ in mesh.vertex(vid).incoming_halfedges() {
            count += 1;
            if count > v_in_count[vid] {
                return Err(format!("vertex {} incomming loop", vid.0));
            }
        }

        if count != v_in_count[vid] {
            return Err(format!("vertex {} incomming loop", vid.0));
        }
    }

    Ok(())
}

#[test]
fn build_nomanifold_mesh() {
    let mesh = SurfaceMesh::new(vec![
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
        .incoming_halfedges()
        .map(|hid| *mesh.he_from(*hid))
        .collect::<Vec<_>>();
    incoming_vertices.sort();
    assert_eq!(base_vertices, incoming_vertices);

    // vertices of outgoing halfedges
    let mut outgoing_vertices = mesh
        .vertex(0.into())
        .outgoing_halfedges()
        .map(|hid| *mesh.he_to(*hid))
        .collect::<Vec<_>>();
    outgoing_vertices.sort();
    assert_eq!(outgoing_vertices, base_vertices);

    // adjacent vertices
    let mut adjacent_vertices = mesh
        .vertex(0.into())
        .vertices()
        .map(|vid| **vid)
        .collect::<Vec<_>>();
    adjacent_vertices.sort();
    assert_eq!(adjacent_vertices, base_vertices);
}

#[test]
fn split_edge_and_face() {
    use bumpalo::collections::Vec;
    let bump = bumpalo::Bump::new();
    let mut mesh = SurfaceMesh::new(vec![
        vec![0, 1, 2],
        vec![0, 1, 3],
        vec![0, 1, 4],
        vec![0, 1, 5],
    ]);
    // split edge
    {
        mesh.split_edge(1.into(), &bump);
        assert!(validate_mesh_connectivity(&mesh).is_ok());
        for f in mesh.faces() {
            let fid = *f;
            let f_verts = Vec::from_iter_in(
                mesh.face(fid).halfedges().map(|hid| mesh.he_from(*hid)),
                &bump,
            );
            assert_eq!(f_verts.len(), 4);
        }

        let new_v = 6.into();
        let mut n_new_halfedges = 0;
        for hid in mesh.vertex(new_v).incoming_halfedges() {
            assert_eq!(mesh.he_to(*hid), new_v);
            n_new_halfedges += 1;
        }
        assert_eq!(n_new_halfedges, 4);
        n_new_halfedges = 0;
        for hid in mesh.vertex(new_v).outgoing_halfedges() {
            assert_eq!(mesh.he_from(*hid), new_v);
            n_new_halfedges += 1;
        }
        assert_eq!(n_new_halfedges, 4);
    }
    let all_edges = Vec::from_iter_in(mesh.edges().map(|e| *e), &bump);
    for eid in all_edges {
        mesh.split_edge(eid, &bump);
        assert!(validate_mesh_connectivity(&mesh).is_ok());
    }
    // split face
    {
        let new_hid = mesh.split_face(0.into(), 2.into(), 6.into(), &bump);
        assert_eq!(
            mesh.face(mesh.he_face(new_hid)).halfedges().count(),
            5,
        );
    }
}
