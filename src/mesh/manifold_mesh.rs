use hashbrown::HashMap;
use itertools::Itertools;

use super::{mesh_core_data::MeshCoreData, EdgeId, ElementId, FaceId, HalfedgeId, Mesh, VertexId};

pub struct ManifoldMesh {
    core_data: MeshCoreData,

    he_face_arr: Vec<FaceId>,
}

impl ManifoldMesh {
    pub fn new<U, T>(polygons: T) -> Self
    where
        U: AsRef<[usize]>,
        T: Iterator<Item = U> + Clone,
    {
        let mut n_faces = 0;
        let mut n_vertices = 0;
        for poly in polygons.clone() {
            n_vertices = n_vertices.max(*poly.as_ref().iter().max().unwrap());
            n_faces += 1;
        }
        n_vertices += 1;
        let core_data = MeshCoreData::new(n_vertices, n_faces);
        let mut mesh = Self {
            core_data,
            he_face_arr: Vec::new(),
        };

        let mut edge_map = HashMap::<(usize, usize), HalfedgeId>::new();
        for (fid, polygon) in polygons.enumerate() {
            let mut first_hid = HalfedgeId::default();
            let mut prev_hid = first_hid;

            for (&a, &b) in polygon.as_ref().iter().circular_tuple_windows() {
                let key = if a < b { (a, b) } else { (b, a) };
                let hid = match edge_map.entry(key) {
                    hashbrown::hash_map::Entry::Occupied(entry) => {
                        let hid = mesh.he_twin(*entry.get());
                        *entry.into_mut() = hid;
                        hid
                    }
                    hashbrown::hash_map::Entry::Vacant(entry) => {
                        let hid = mesh.new_edge();
                        entry.insert(hid);
                        hid
                    }
                };
                mesh.core_data.v_halfedge_arr[a] = hid;
                mesh.core_data.set_he_vertex(hid, b.into());
                mesh.core_data.set_he_vertex(mesh.he_twin(hid), a.into());
                mesh.he_face_arr[hid] = fid.into();
                if !first_hid.valid() {
                    mesh.core_data.f_halfedge_arr[fid] = hid;
                    first_hid = hid;
                } else {
                    mesh.core_data.connect_halfedges(prev_hid, hid);
                }
                prev_hid = hid;
            }
            mesh.core_data.connect_halfedges(prev_hid, first_hid);
        }

        let mut edge_visited = vec![false; mesh.n_edges_capacity()];
        for twin_hid in edge_map.into_values() {
            if twin_hid.0 & 1 != 0 {
                continue;
            }
            let eid = mesh.he_edge(twin_hid);
            if edge_visited[eid] {
                continue;
            }
            edge_visited[eid] = true;

            let first_hid = mesh.he_twin(twin_hid);
            debug_assert!(mesh.he_is_boundary(first_hid));
            let mut curr_hid = first_hid;
            loop {
                let mut prev_hid = curr_hid;
                loop {
                    prev_hid = mesh.he_next(mesh.he_twin(prev_hid));
                    let twin_prev_hid = mesh.he_twin(prev_hid);
                    if mesh.he_is_boundary(twin_prev_hid) {
                        prev_hid = twin_prev_hid;
                        break;
                    }
                }
                mesh.core_data
                    .set_v_halfedge(mesh.he_to(prev_hid), curr_hid);
                edge_visited[mesh.he_edge(prev_hid)] = true;
                mesh.core_data.connect_halfedges(prev_hid, curr_hid);
                curr_hid = prev_hid;
                if curr_hid == first_hid {
                    break;
                }
            }
        }
        mesh
    }

    pub fn new_face_by_halfedges(&mut self, halfedges: &[HalfedgeId]) -> FaceId {
        for (&ha_twin, &hb_twin) in halfedges.iter().rev().circular_tuple_windows() {
            let ha = self.he_twin(ha_twin);
            let hb = self.he_twin(hb_twin);
            match [self.he_is_boundary(ha), self.he_is_boundary(hb)] {
                [true, true] => {
                    self.core_data.connect_halfedges(ha, hb);
                    self.core_data.set_v_halfedge(self.he_to(ha), hb);
                }
                [true, false] => {
                    let ha_next = self.he_next(hb_twin);
                    self.core_data.connect_halfedges(ha, ha_next);
                }
                [false, true] => {
                    let hb_prev = self.he_prev(ha_twin);
                    self.core_data.connect_halfedges(hb_prev, hb);
                    self.core_data.set_v_halfedge(self.he_to(ha), hb);
                }
                [false, false] => {}
            }
        }
        let new_fid = self.new_face();
        for (&ha, &hb) in halfedges.iter().circular_tuple_windows() {
            self.core_data.connect_halfedges(ha, hb);
            self.he_face_arr[ha] = new_fid;
        }
        self.core_data.f_halfedge_arr[new_fid] = halfedges[0];
        new_fid
    }

    #[inline(always)]
    pub fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.he_face_arr[hid]
    }

    #[inline(always)]
    pub fn he_is_boundary(&self, hid: HalfedgeId) -> bool {
        !self.he_face(hid).valid()
    }

    #[inline(always)]
    pub fn new_vertices(&mut self, len: usize) -> VertexId {
        let new_vid = VertexId::from(self.n_vertices_capacity());
        let new_len = new_vid.0 + len;
        self.core_data
            .v_halfedge_arr
            .resize(new_len, HalfedgeId::default());
        self.core_data.n_vertices += len;
        new_vid
    }

    #[inline]
    pub fn new_edge(&mut self) -> HalfedgeId {
        let hid = self.core_data.he_next_arr.len().into();
        let new_len = self.core_data.he_prev_arr.len() + 2;

        self.core_data
            .he_prev_arr
            .resize(new_len, HalfedgeId::default());
        self.core_data
            .he_next_arr
            .resize(new_len, HalfedgeId::default());
        self.core_data
            .he_vertex_arr
            .resize(new_len, VertexId::default());

        self.he_face_arr.resize(new_len, FaceId::default());
        self.core_data.n_halfedges += 2;
        hid
    }

    #[inline]
    pub fn new_face(&mut self) -> FaceId {
        let fid = FaceId::from(self.core_data.f_halfedge_arr.len());
        self.core_data.f_halfedge_arr.push(HalfedgeId::default());
        self.core_data.n_faces += 1;
        fid
    }

    #[inline]
    pub fn new_edge_by_veritces(&mut self, va: VertexId, vb: VertexId) -> HalfedgeId {
        let hid = self.new_edge();
        let twin_hid = self.he_twin(hid);
        self.core_data.set_he_vertex(hid, vb);
        self.core_data.set_he_vertex(twin_hid, va);
        hid
    }

    #[inline]
    pub fn remove_vertex(&mut self, vid: VertexId) {
        self.core_data.set_v_halfedge(vid, HalfedgeId::default());
        self.core_data.n_vertices -= 1;
    }

    #[inline]
    pub fn remove_edge(&mut self, eid: EdgeId) {
        let hid = self.e_halfedge(eid);
        self.core_data.set_he_vertex(hid, VertexId::default());
        self.core_data
            .set_he_vertex(self.he_twin(hid), VertexId::default());
        self.core_data.n_halfedges -= 2;
    }

    pub fn remove_face(&mut self, fid: FaceId) {
        let first_hid = self.f_halfedge(fid);
        let prev_first_hid = self.he_prev(first_hid);
        let mut curr_hid = first_hid;
        loop {
            self.he_face_arr[curr_hid] = FaceId::default();

            let next_hid = if curr_hid == prev_first_hid {
                first_hid
            } else {
                self.he_next(curr_hid)
            };
            let rev_next_hid = self.he_twin(curr_hid);
            let rev_curr_hid = self.he_twin(next_hid);

            match [
                self.he_is_boundary(rev_next_hid),
                self.he_is_boundary(rev_curr_hid),
            ] {
                [true, true] => {
                    // \  /
                    //  \/
                    //  vid
                    //  /\
                    // /  \
                    let vid = self.he_to(curr_hid);
                    let vh = self.v_halfedge(vid);

                    if vh == rev_next_hid && self.he_prev(rev_next_hid) == rev_curr_hid {
                        self.remove_vertex(vid);
                    } else {
                        let prev_rev_next_hid = self.he_prev(rev_next_hid);
                        debug_assert!(self.he_is_valid(prev_rev_next_hid));
                        let next_rev_curr_hid = self.he_next(rev_curr_hid);
                        debug_assert!(self.he_is_valid(next_rev_curr_hid));
                        self.core_data
                            .connect_halfedges(prev_rev_next_hid, next_rev_curr_hid);
                        self.core_data.set_v_halfedge(vid, next_rev_curr_hid);
                    }
                }
                [true, false] => {
                    self.core_data
                        .connect_halfedges(self.he_prev(rev_next_hid), next_hid);
                    self.core_data
                        .set_v_halfedge(self.he_to(curr_hid), next_hid);
                }
                [false, true] => {
                    let next_rev_curr_hid = self.he_next(rev_curr_hid);
                    self.core_data
                        .connect_halfedges(curr_hid, next_rev_curr_hid);
                    self.core_data
                        .set_v_halfedge(self.he_to(curr_hid), next_rev_curr_hid);
                }
                [false, false] => {
                    self.core_data
                        .set_v_halfedge(self.he_to(curr_hid), next_hid);
                }
            }

            // make current edge invalid
            if self.he_is_boundary(rev_next_hid) {
                self.core_data.set_he_vertex(curr_hid, VertexId::default());
                self.core_data
                    .set_he_vertex(rev_next_hid, VertexId::default());
                self.core_data.n_halfedges -= 2;
            }

            curr_hid = next_hid;
            if curr_hid == first_hid {
                break;
            }
        }
        self.core_data.f_halfedge_arr[fid] = HalfedgeId::default();
        self.core_data.n_faces -= 1;
    }
}

impl Mesh for ManifoldMesh {
    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        true
    }

    #[inline(always)]
    fn core_data(&self) -> &MeshCoreData {
        &self.core_data
    }

    #[inline(always)]
    fn n_edges(&self) -> usize {
        self.core_data.n_halfedges >> 1
    }

    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        self.n_halfedges_capacity() >> 1
    }

    #[inline(always)]
    fn e_is_valid(&self, eid: EdgeId) -> bool {
        let idx = eid.0 << 1;
        self.he_is_valid(idx.into()) || self.he_is_valid((idx + 1).into())
    }

    #[inline(always)]
    fn he_from(&self, hid: HalfedgeId) -> VertexId {
        self.he_to(self.he_twin(hid))
    }

    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        (hid.0 ^ 1).into()
    }

    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_twin(hid)
    }

    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_prev(self.he_twin(hid))
    }

    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        (hid.0 >> 1).into()
    }

    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        (eid.0 << 1).into()
    }
}

pub fn validate_mesh_connectivity(mesh: &ManifoldMesh) -> Result<(), String> {
    let validate_vertex = |vid: VertexId, msg: &str| {
        if vid.0 > mesh.n_vertices_capacity() || !mesh.v_is_valid(vid) {
            Err(format!("{} bad vertex reference: {}", msg, vid.0))
        } else {
            Ok(())
        }
    };
    let validate_halfedge = |hid: HalfedgeId, msg: &str| {
        if hid.0 > mesh.n_halfedges_capacity() || !mesh.he_is_valid(hid) {
            Err(format!("{} bad halfedge reference: {}", msg, hid.0))
        } else {
            Ok(())
        }
    };
    let validate_edge = |eid: EdgeId, msg: &str| {
        if eid.0 > mesh.n_edges_capacity() || !mesh.e_is_valid(eid) {
            Err(format!("{} bad edge reference: {}", msg, eid.0))
        } else {
            Ok(())
        }
    };
    let validate_face = |fid: FaceId, msg: &str| {
        if fid.0 > mesh.n_faces_capacity() || !mesh.f_is_valid(fid) {
            Err(format!("{} bad face reference: {}", msg, fid.0))
        } else {
            Ok(())
        }
    };

    for v in mesh.vertices() {
        validate_halfedge(mesh.v_halfedge(*v), "v_he: ")?;
    }

    for he in mesh.halfedges() {
        validate_vertex(mesh.he_from(*he), "he_vertex: ")?;

        validate_halfedge(mesh.he_next(*he), "he_next: ")?;
        validate_halfedge(mesh.he_twin(*he), "he_twin: ")?;
        validate_halfedge(mesh.he_next_incoming_neighbor(*he), "next_incoming: ")?;

        validate_edge(mesh.he_edge(*he), "he_edge: ")?;
    }

    for e in mesh.edges() {
        validate_halfedge(mesh.e_halfedge(*e), "e_he: ")?;
    }

    for f in mesh.faces() {
        validate_halfedge(mesh.f_halfedge(*f), "e_face: ")?;
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

    for he in mesh.halfedges() {
        let hid = *he;
        let tip = mesh.he_to(hid);
        if *mesh.halfedge(mesh.he_next_incoming_neighbor(hid)).to() != tip {
            return Err(format!("next incoming he is not to same vert"));
        }
    }

    for he in mesh.halfedges() {
        let next_he = he.next();
        if *he.to() != *next_he.from() {
            return Err(format!(
                "prev he {} and next he {} are not connected",
                he.0, next_he.0
            ));
        }
    }

    for he in mesh.halfedges() {
        if mesh.he_is_boundary(*he) {
            if !mesh.he_is_boundary(mesh.v_halfedge(*he.from())) {
                return Err(format!(
                    "the halfedge of  boundary vertex {:?} is not boundary",
                    *he.from()
                ));
            }
        }
    }
    {
        let mut visisted = vec![false; mesh.n_halfedges_capacity()];
        for he in mesh.halfedges() {
            if visisted[*he] {
                continue;
            }

            let first = *he;
            let mut curr = first;
            let mut count = 0;
            {
                count += 1;
                if count > mesh.n_halfedges_capacity() {
                    return Err(format!("the loop from {:?} is not closed", first));
                }
                visisted[curr] = true;
                curr = mesh.he_next(curr);
                if curr == first {
                    break;
                }
            }
        }
    }

    let mut v_in_count = vec![0; mesh.n_vertices_capacity()];
    for hid in mesh.halfedges() {
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
