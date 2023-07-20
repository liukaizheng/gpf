mod element;
mod manifold_mesh;
mod mesh_data;
mod surface_mesh;

#[macro_use]
mod mesh_macro;

use std::{cell::RefCell, rc::Weak};

pub use element::*;
pub use manifold_mesh::*;
pub use mesh_data::*;
pub use surface_mesh::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FaceOrBoundaryLoopId {
    Face(FaceId),
    BoundaryLoop(BoundaryLoopId),
    INVALID,
}

pub trait Mesh: Sized {
    fn n_vertices(&self) -> usize;
    fn n_halfedges(&self) -> usize;
    fn n_edges(&self) -> usize;
    fn n_faces(&self) -> usize;
    // fn n_boundary_loops(&self) -> usize;

    fn n_vertices_capacity(&self) -> usize;
    fn n_halfedges_capacity(&self) -> usize;
    fn n_edges_capacity(&self) -> usize;
    fn n_faces_capacity(&self) -> usize;
    // fn n_boundary_loop_capacity(&self) -> usize;

    fn vertex_is_valid(&self, vid: VertexId) -> bool;
    fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool;
    fn edge_is_valid(&self, eid: EdgeId) -> bool;
    fn face_is_valid(&self, fid: FaceId) -> bool;
    // fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool;

    fn vertex(&self, vid: VertexId) -> VertexIter<Self>;
    fn halfedge(&self, hid: HalfedgeId) -> HalfedgeIter<Self>;
    fn edge(&self, eid: EdgeId) -> EdgeIter<Self>;
    fn face(&self, fid: FaceId) -> FaceIter<Self>;
    // fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self>;

    fn vertices(&self) -> VertexIter<Self>;
    fn halfedges(&self) -> HalfedgeIter<Self>;
    fn edges(&self) -> EdgeIter<Self>;
    fn faces(&self) -> FaceIter<Self>;

    /// the halfedge starting from this vertex
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId;

    /// the start vertex of the halfedge
    fn he_vertex(&self, hid: HalfedgeId) -> VertexId;
    /// the end vertex of the halfedge
    fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId;
    /// the next halfedge of the halfedge
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the previous halfedge of the halfedge
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the twin halfedge of the halfedge
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the sibling halfedge of the halfedge
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the next incoming halfedge of the halfedge
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the next outgoing halfedge of the halfedge
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId;
    /// the edge of the halfedge
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId;
    /// the face or boundary loop of the halfedge
    fn he_face_or_boundary_loop(&self, hid: HalfedgeId) -> FaceOrBoundaryLoopId;

    /// the first halfedge of the edge
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId;

    /// two vertices of the edge
    fn e_vertices(&self, eid: EdgeId) -> [VertexId; 2];

    /// the first halfedge of the face
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId;

    #[inline(always)]
    fn eh_same_dir(&self, eid: EdgeId, hid: HalfedgeId) -> bool {
        self.he_vertex(hid) == self.he_vertex(self.e_halfedge(eid))
    }

    fn use_implicit_twin(&self) -> bool;
    /// add vertices data
    fn add_vertices_data<T: Clone + 'static>(&mut self, data: Weak<RefCell<VertexData<T, Self>>>);
    fn remove_vertices_data<T: Clone + 'static>(&mut self, data: &VertexData<T, Self>);

    /// add halfedges data reference
    fn add_halfedges_data<T: Clone + 'static>(
        &mut self,
        data: Weak<RefCell<HalfedgeData<T, Self>>>,
    );
    fn remove_halfedges_data<T: Clone + 'static>(&mut self, data: &HalfedgeData<T, Self>);

    fn add_edges_data<T: Clone + 'static>(&mut self, data: Weak<RefCell<EdgeData<T, Self>>>);
    fn remove_edges_data<T: Clone + 'static>(&mut self, data: &EdgeData<T, Self>);

    fn add_faces_data<T: Clone + 'static>(&mut self, data: Weak<RefCell<FaceData<T, Self>>>);
    fn remove_faces_data<T: Clone + 'static>(&mut self, data: &FaceData<T, Self>);
}

pub fn validate_mesh_connectivity(mesh: &SurfaceMesh) -> Result<(), String> {
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
        validate_halfedge(mesh.v_halfedge(vid), "v_he: ")?;
    }

    for hid in mesh.halfedges() {
        validate_vertex(mesh.he_vertex(hid), "he_vertex: ")?;

        validate_halfedge(mesh.he_next(hid), "he_next: ")?;
        validate_halfedge(mesh.he_twin(hid), "he_twin: ")?;
        validate_halfedge(
            mesh.he_next_incoming_neighbor(hid),
            &format!("next_incoming for {}: ", hid.0),
        )?;
        validate_halfedge(mesh.he_next_outgoing_neighbor(hid), "next_outgoing: ")?;

        validate_edge(mesh.he_edge(hid), "he_edge: ")?;

        if let FaceOrBoundaryLoopId::Face(fid) = mesh.he_face_or_boundary_loop(hid) {
            validate_face(fid, "he_face: ")?;
        }
    }

    for eid in mesh.edges() {
        validate_halfedge(mesh.e_halfedge(eid), "e_he: ")?;
    }

    for fid in mesh.faces() {
        validate_halfedge(mesh.f_halfedge(fid), "e_face: ")?;
    }

    for hid in mesh.halfedges() {
        let first_he = mesh.halfedge(hid);
        let mut curr_he = mesh.halfedge(hid);
        curr_he.sibling();
        let mut count = 1;
        while *curr_he != hid {
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
            curr_he.sibling();
            count += 1;
        }
    }

    for eid in mesh.edges() {
        for hid in mesh.edge(eid).halfedges() {
            if eid != mesh.he_edge(hid) {
                return Err(format!(
                    "edge.halfedge doesn't match he.edge for edge {}",
                    eid.0
                ));
            }
        }
    }

    for fid in mesh.faces() {
        let hid = mesh.f_halfedge(fid);
        if mesh.he_face_arr[hid] != fid {
            return Err(format!("f.he().face() is not f for face {}", fid.0));
        }

        let mut curr_he = mesh.halfedge(hid);
        Halfedge::next(&mut curr_he);
        let mut count = 1;
        while *curr_he != hid {
            if *curr_he.face().unwrap() != fid {
                return Err(format!("face.he doesn't match he.face for face {}", fid.0));
            }
            if count > mesh.n_halfedges() {
                return Err(format!(
                    "halfedge next doesn't cycle back for face {}",
                    fid.0
                ));
            }
            Halfedge::next(&mut curr_he);
            count += 1;
        }
        if count < 2 {
            return Err(format!("face {} degree < 2", fid.0));
        }
    }

    for hid in mesh.halfedges() {
        let tail = mesh.he_vertex(hid);
        let tip = mesh.he_tip_vertex(hid);
        if *mesh
            .halfedge(mesh.he_next_incoming_neighbor(hid))
            .tip_vertex()
            != tip
        {
            return Err(format!("next incoming he is not to same vert"));
        }
        if *mesh.halfedge(mesh.he_next_outgoing_neighbor(hid)).vertex() != tail {
            return Err(format!("next outgoing he is not from same vert"));
        }
    }

    let bump = bumpalo::Bump::new();
    let mut v_in_count = bumpalo::vec![in &bump; 0; mesh.n_vertices_capacity()];
    let mut v_out_count = bumpalo::vec![in &bump; 0; mesh.n_vertices_capacity()];
    for hid in mesh.halfedges() {
        v_out_count[mesh.he_vertex(hid)] += 1;
        v_in_count[mesh.he_tip_vertex(hid)] += 1;
    }
    for vid in mesh.vertices() {
        let mut count = 0;
        for _ in mesh.vertex(vid).incoming_halfedge() {
            count += 1;
            if count > v_in_count[vid] {
                return Err(format!("vertex {} incomming loop", vid.0));
            }
        }

        if count != v_in_count[vid] {
            return Err(format!("vertex {} incomming loop", vid.0));
        }

        count = 0;
        for _ in mesh.vertex(vid).outgoing_halfedge() {
            count += 1;
            if count > v_out_count[vid] {
                return Err(format!("vertex {} outgoing loop", vid.0));
            }
        }

        if count != v_out_count[vid] {
            return Err(format!("vertex {} outgoing loop", vid.0));
        }
    }

    Ok(())
}
