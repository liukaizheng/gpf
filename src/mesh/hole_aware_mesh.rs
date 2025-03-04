use std::alloc::Allocator;

use hashbrown::HashMap;

use super::{mesh_core_data::MeshCoreData, EdgeId, ElementId, FaceId, HalfedgeId, Mesh, VertexId};

pub struct HoleAwareMesh<A: Allocator + Copy> {
    core_data: MeshCoreData<A>,
    n_edges: usize,
    he_edge_arr: Vec<EdgeId, A>,
    he_vert_in_next_arr: Vec<HalfedgeId, A>,
    he_face_next_loop_arr: Vec<HalfedgeId, A>,
    he_sibling_arr: Vec<HalfedgeId, A>,
    e_halfedge_arr: Vec<HalfedgeId, A>,
}

impl<A: Allocator + Copy> HoleAwareMesh<A> {
    pub fn new<U1, T1, U2, T2>(loops: T1, faces: T2, alloc: A) -> Self
    where
        U1: AsRef<[usize]>,
        T1: IntoIterator<Item = U1>,
        U2: AsRef<[usize]>,
        T2: IntoIterator<Item = U2>,
    {
        let core_data = MeshCoreData::new(0, 0, alloc);
        let mut mesh = Self {
            core_data,
            n_edges: 0,
            he_edge_arr: Vec::new_in(alloc),
            he_vert_in_next_arr: Vec::new_in(alloc),
            he_face_next_loop_arr: Vec::new_in(alloc),
            he_sibling_arr: Vec::new_in(alloc),
            e_halfedge_arr: Vec::new_in(alloc),
        };

        let mut loop_halfedges: Vec<HalfedgeId, A> = Vec::new_in(alloc);
        for loop_vertices in loops.into_iter() {
            let mut first_hid = HalfedgeId::default();
            let mut prev_vid = VertexId::default();
            let mut prev_hid = HalfedgeId::default();
            for (i, &vid) in loop_vertices.as_ref().iter().enumerate() {
                let vid = VertexId(vid);
                mesh.core_data.v_min_reserve(vid);
                let hid = mesh.new_halfedges(1);

                mesh.core_data.he_vertex_arr[hid] = vid;

                if i == 0 {
                    first_hid = hid;
                } else {
                    mesh.core_data.v_halfedge_arr[prev_vid] = hid;
                    mesh.core_data.connect_halfedges(prev_hid, hid);
                }
                prev_vid = vid;
                prev_hid = hid;
            }
            mesh.core_data.v_halfedge_arr[prev_vid] = first_hid;
            mesh.core_data.connect_halfedges(prev_hid, first_hid);
            loop_halfedges.push(first_hid);
        }

        for (fid, loop_indices) in faces.into_iter().enumerate() {
            let fid = FaceId(fid);
            let mut first_loop_hid = HalfedgeId::default();
            let mut prev_loop_hid = HalfedgeId::default();
            for (i, &lid) in loop_indices.as_ref().iter().enumerate() {
                let curr_loop_hid = loop_halfedges[lid];
                if i == 0 {
                    first_loop_hid = curr_loop_hid;
                    prev_loop_hid = first_loop_hid;
                } else {
                    mesh.he_face_next_loop_arr[prev_loop_hid] = curr_loop_hid;
                    prev_loop_hid = curr_loop_hid;
                }

                // set face id for all halfedges in the loop
                let mut curr_hid = curr_loop_hid;
                loop {
                    mesh.core_data.set_he_face(curr_hid, fid);
                    curr_hid = mesh.he_next(curr_hid);
                    if curr_hid == curr_loop_hid {
                        break;
                    }
                }
            }
            mesh.he_face_next_loop_arr[prev_loop_hid] = first_loop_hid;
            mesh.core_data.f_halfedge_arr.push(first_loop_hid);
        }
        mesh.core_data.n_faces = mesh.core_data.f_halfedge_arr.len();
        mesh.core_data.recount_n_vertices();

        let mut edge_history = HashMap::<(usize, usize), HalfedgeId>::new();
        // build edge
        for hid in 0..mesh.he_edge_arr.len() {
            let hid = hid.into();
            let [va, vb] = mesh.he_vertices(hid);
            let key = if va.0 < vb.0 {
                (va.0, vb.0)
            } else {
                (vb.0, va.0)
            };
            if let Some(prev_hid) = edge_history.get_mut(&key) {
                // We're already seen this edge, connect to the previous halfedge incident on the edge
                mesh.he_sibling_arr[hid] = *prev_hid;
                let eid = mesh.he_edge_arr[*prev_hid];
                mesh.he_edge_arr[hid] = eid;
                *prev_hid = hid;
            } else {
                // This is the first time we've ever seen this edge, create a new edge object
                let new_eid = mesh.new_edges(1);
                mesh.he_edge_arr[hid] = new_eid;
                mesh.he_sibling_arr[hid] = HalfedgeId::default();
                mesh.e_halfedge_arr[new_eid] = hid;
                edge_history.insert(key, hid);
            }
        }
        // Complete the sibling cycle by following backwards each edge until we reach the first sibling-less entry
        for last_he in edge_history.into_values() {
            if !mesh.he_sibling_arr[last_he].valid() {
                // Any edges which never got any sibling entries at all are boundary halfedges
                mesh.he_sibling_arr[last_he] = last_he;
                continue;
            }

            // Get the index of the first halfedge in the sibling cycle to complete the cycle
            let mut curr_he = last_he;
            while mesh.he_sibling_arr[curr_he].valid() {
                curr_he = mesh.he_sibling_arr[curr_he];
            }
            mesh.he_sibling_arr[curr_he] = last_he; // connect the first to the last
        }

        let (v_in_halfedges, v_in_separators) = mesh.vertex_cycle();
        let n_vertices = mesh.n_vertices();
        for idx in 0..n_vertices {
            let vid = VertexId::from(idx);
            if !mesh.v_is_valid(vid) {
                continue;
            }
            let (start, end) = (v_in_separators[vid], v_in_separators[*vid + 1]);
            let len = end - start;
            for i in start..end {
                let ha = v_in_halfedges[i];
                let hb = v_in_halfedges[start + (i - start + 1) % len];
                mesh.he_vert_in_next_arr[ha] = hb;
            }
        }

        mesh
    }

    fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        let hid = HalfedgeId::from(self.core_data.he_next_arr.len());
        let new_len = self.n_halfedges_capacity() + n;
        self.core_data.new_halfedges(n);
        self.he_sibling_arr.resize(new_len, HalfedgeId::default());
        self.he_edge_arr.resize(new_len, EdgeId::default());
        self.he_vert_in_next_arr
            .resize(new_len, HalfedgeId::default());
        self.he_face_next_loop_arr
            .resize(new_len, HalfedgeId::default());

        self.core_data.n_halfedges += n;
        hid
    }

    #[inline]
    pub fn new_edges(&mut self, n: usize) -> EdgeId {
        let eid = self.e_halfedge_arr.len().into();
        let new_len = self.n_edges_capacity() + n;
        self.e_halfedge_arr.resize(new_len, HalfedgeId::default());

        self.n_edges += n;
        eid
    }

    #[inline]
    fn vertex_cycle(&self) -> (Vec<HalfedgeId>, Vec<usize>) {
        let mut v_degree = vec![0usize; self.n_vertices_capacity()];
        self.halfedges().for_each(|he| {
            let vertex = he.to();
            v_degree[*vertex] += 1;
        });
        let mut vertex_separators = vec![0];
        vertex_separators.extend(v_degree.iter().scan(0, |sum, &count| {
            *sum += count;
            Some(*sum)
        }));
        let mut he_positions = vertex_separators.clone();
        let mut vertex_halfedges = vec![HalfedgeId::from(0); self.n_halfedges_capacity()];
        self.halfedges().for_each(|he| {
            let vid = *he.to();
            let pos = he_positions[vid];
            vertex_halfedges[pos] = *he;
            he_positions[vid] += 1;
        });
        (vertex_halfedges, vertex_separators)
    }
}

impl<A: Allocator + Copy> Mesh for HoleAwareMesh<A> {
    #[inline]
    fn n_vertices(&self) -> usize {
        self.core_data.n_vertices
    }

    #[inline]
    fn n_halfedges(&self) -> usize {
        self.core_data.n_halfedges
    }

    #[inline]
    fn n_edges(&self) -> usize {
        self.n_edges
    }

    #[inline]
    fn n_faces(&self) -> usize {
        self.core_data.n_faces
    }

    #[inline]
    fn n_vertices_capacity(&self) -> usize {
        self.core_data.n_vertices_capacity()
    }

    #[inline]
    fn n_halfedges_capacity(&self) -> usize {
        self.core_data.n_halfedges_capacity()
    }

    #[inline]
    fn n_edges_capacity(&self) -> usize {
        self.e_halfedge_arr.len()
    }

    #[inline]
    fn n_faces_capacity(&self) -> usize {
        self.core_data.n_faces_capacity()
    }

    fn set_n_vertices(&mut self, n: usize) {
        self.core_data.n_vertices = n;
        if n > 0 {
            self.core_data.v_min_reserve((n - 1).into());
        }
    }

    #[inline]
    fn e_is_valid(&self, eid: EdgeId) -> bool {
        self.e_halfedge(eid).valid()
    }

    #[inline(always)]
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
        self.core_data.v_halfedge_arr[vid]
    }

    #[inline(always)]
    fn he_to(&self, hid: HalfedgeId) -> VertexId {
        self.core_data.he_vertex_arr[hid]
    }

    #[inline]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.core_data.he_next_arr[hid]
    }

    #[inline]
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
        self.core_data.he_prev_arr[hid]
    }

    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        let mut curr = hid;
        loop {
            let next = self.he_sibling_arr[curr];
            if next == hid {
                return HalfedgeId::default();
            }

            if self.he_to(curr) != self.he_to(next) {
                return next;
            }

            curr = next;
        }
    }

    #[inline]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_sibling_arr[hid]
    }

    #[inline]
    fn he_next_incoming(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_vert_in_next_arr[hid]
    }

    #[inline]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        self.he_edge_arr[hid]
    }

    #[inline]
    fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.core_data.he_face_arr[hid]
    }

    #[inline]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        self.e_halfedge_arr[eid]
    }

    #[inline]
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
        self.core_data.f_halfedge_arr[fid]
    }

    #[inline(always)]
    fn f_loop_next_first_halfedge(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_face_next_loop_arr[hid]
    }

    #[inline]
    fn use_implicit_twin(&self) -> bool {
        false
    }

    #[inline]
    fn set_v_halfedge(&mut self, v: VertexId, hid: HalfedgeId) {
        self.core_data.set_v_halfedge(v, hid);
    }

    #[inline]
    fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId) {
        self.core_data.set_f_halfedge(fid, hid);
    }

    #[inline]
    fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId) {
        self.core_data.set_he_vertex(hid, vid);
    }
}
