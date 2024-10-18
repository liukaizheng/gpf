use std::alloc::Allocator;

use crate::mesh::Mesh;

use super::element::{EdgeId, ElementId, FaceId, HalfedgeId, VertexId};
use super::mesh_core_data::MeshCoreData;
use hashbrown::HashMap;
use itertools::Itertools;

pub struct SurfaceMesh {
    core_data: MeshCoreData,
    n_edges: usize,
    he_edge_arr: Vec<EdgeId>,
    he_face_arr: Vec<FaceId>,
    he_vert_in_next_arr: Vec<HalfedgeId>,
    he_sibling_arr: Vec<HalfedgeId>,
    e_halfedge_arr: Vec<HalfedgeId>,

    use_implicit_twin: bool,
}

impl SurfaceMesh {
    pub fn new<T, U>(polygons: T) -> Self
    where
        T: AsRef<[U]>,
        U: AsRef<[usize]>,
    {
        let n_faces = polygons.as_ref().len();
        let n_vertices = polygons
            .as_ref()
            .iter()
            .map(|polygon| *polygon.as_ref().iter().max().unwrap())
            .max()
            .unwrap()
            + 1;
        let core_data = MeshCoreData::new(n_vertices, n_faces);
        let mut mesh = Self {
            core_data,
            n_edges: 0,
            he_edge_arr: Vec::new(),
            he_face_arr: Vec::new(),
            he_vert_in_next_arr: Vec::new(),
            he_sibling_arr: Vec::new(),
            e_halfedge_arr: Vec::new(),
            use_implicit_twin: false,
        };

        for (fid, polygon) in polygons.as_ref().iter().enumerate() {
            let mut first_hid = HalfedgeId::default();
            let mut prev_hid = first_hid;
            let mut prev_vid = VertexId::default();

            for (i, &b) in polygon.as_ref().iter().enumerate() {
                let vid = b.into();
                let hid = mesh.new_halfedges(1);

                mesh.core_data.he_vertex_arr[hid] = vid;
                mesh.he_face_arr[hid] = fid.into();

                if i == 0 {
                    mesh.core_data.f_halfedge_arr[fid] = hid;
                    first_hid = hid;
                } else {
                    mesh.core_data.v_halfedge_arr[prev_vid] = hid;
                    mesh.core_data.he_next_arr[prev_hid] = hid;
                    mesh.core_data.connect_halfedges(prev_hid, hid);
                }
                prev_vid = vid;
                prev_hid = hid;
            }
            mesh.core_data.v_halfedge_arr[prev_vid] = first_hid;
            mesh.core_data.connect_halfedges(prev_hid, first_hid);
        }

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
            if !mesh.vertex_is_valid(vid) {
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

    #[inline]
    pub fn n_edges_capacity(&self) -> usize {
        self.e_halfedge_arr.len()
    }

    #[inline]
    pub fn new_vertices(&mut self, n: usize) -> VertexId {
        let vid = VertexId::from(self.core_data.v_halfedge_arr.len());
        let new_len = self.n_vertices_capacity() + n;
        self.core_data
            .v_halfedge_arr
            .resize(new_len, HalfedgeId::default());
        self.core_data.n_vertices += n;
        vid
    }

    pub fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        let hid = HalfedgeId::from(self.core_data.he_next_arr.len());
        let new_len = self.n_halfedges_capacity() + n;
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
        self.he_sibling_arr.resize(new_len, HalfedgeId::default());
        self.he_edge_arr.resize(new_len, EdgeId::default());
        self.he_vert_in_next_arr
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
    pub fn new_faces(&mut self, n: usize) -> FaceId {
        let fid = FaceId::from(self.core_data.f_halfedge_arr.len());
        let new_len = self.core_data.f_halfedge_arr.len() + n;
        self.core_data
            .f_halfedge_arr
            .resize(new_len, HalfedgeId::default());
        self.core_data.n_faces += n;
        fid
    }

    pub fn split_edge<A: Allocator + Copy>(&mut self, eid: EdgeId, bump: A) -> VertexId {
        let mut e_halfedges = Vec::new_in(bump);
        e_halfedges.extend(self.edge(eid).halfedges().map(|he| *he));
        let vb = self.he_to(e_halfedges[0]);
        let new_v = self.new_vertices(1);
        let first_new_he = self.new_halfedges(e_halfedges.len());
        let new_e = self.new_edges(1);

        self.core_data.v_halfedge_arr[new_v] = e_halfedges[0];

        let mut old_e_halfedges = Vec::with_capacity_in(e_halfedges.len(), bump);
        let mut new_e_halfedges = Vec::with_capacity_in(e_halfedges.len(), bump);
        for (&old_hid, idx) in e_halfedges
            .iter()
            .zip(first_new_he.0..self.n_halfedges_capacity())
        {
            let new_hid = idx.into();
            self.core_data.he_vertex_arr[new_hid] = new_v;
            let vid = self.he_from(old_hid);
            if self.v_halfedge(vid) == old_hid {
                self.core_data.v_halfedge_arr[vid] = new_hid;
            }
            if vid != vb {
                self.he_edge_arr[new_hid] = new_e;
                old_e_halfedges.push(old_hid);
                new_e_halfedges.push(new_hid);
            } else {
                self.he_edge_arr[old_hid] = new_e;
                self.he_edge_arr[new_hid] = eid;
                old_e_halfedges.push(new_hid);
                new_e_halfedges.push(old_hid);
            }

            self.he_face_arr[new_hid] = self.he_face_arr[old_hid];

            let prev_hid = self.he_prev(old_hid);
            self.core_data.connect_halfedges(prev_hid, new_hid);
            self.core_data.connect_halfedges(new_hid, old_hid);
        }

        for (ha, hb) in (first_new_he.0..self.n_halfedges_capacity()).circular_tuple_windows() {
            self.he_vert_in_next_arr[ha] = hb.into();
        }

        for (ha, hb) in old_e_halfedges.into_iter().circular_tuple_windows() {
            self.he_sibling_arr[ha] = hb;
        }

        for (ha, hb) in new_e_halfedges.into_iter().circular_tuple_windows() {
            self.he_sibling_arr[ha] = hb;
        }
        new_v
    }

    pub fn split_face<A: Allocator + Copy>(
        &mut self,
        fid: FaceId,
        va: VertexId,
        vb: VertexId,
        bump: A,
    ) -> HalfedgeId {
        let mut right_first_hid = HalfedgeId::default();
        let mut left_first_hid = HalfedgeId::default();
        for he in self.face(fid).halfedges() {
            let v = *he.to();
            if v == va || v == vb {
                if !right_first_hid.valid() {
                    right_first_hid = *he.next();
                } else {
                    left_first_hid = *he.next();
                    break;
                }
            }
        }
        debug_assert!(right_first_hid.valid());
        debug_assert!(left_first_hid.valid());
        let (va, vb) = if self.he_from(right_first_hid) == va {
            (va, vb)
        } else {
            (vb, va)
        };

        //  --------------vb------------
        //  |             |            |
        //  |             |            |
        //  --------------va------------
        let left_last_hid = self.he_prev(right_first_hid);
        let right_last_hid = self.he_prev(left_first_hid);

        let first_he = self.new_halfedges(2);
        let second_he = HalfedgeId::from(first_he.0 + 1);

        // h-v
        self.core_data.he_vertex_arr[first_he] = va;
        self.core_data.he_vertex_arr[second_he] = vb;

        // h-h
        // sibling halfedges
        self.he_sibling_arr[first_he] = second_he;
        self.he_sibling_arr[second_he] = first_he;
        // incoming halfedges and outcoming halfedges
        insert_halfedge(&mut self.he_vert_in_next_arr, left_last_hid, first_he);
        insert_halfedge(&mut self.he_vert_in_next_arr, right_last_hid, second_he);
        // next halfedges
        self.core_data.connect_halfedges(right_last_hid, first_he);
        self.core_data.connect_halfedges(first_he, right_first_hid);
        self.core_data.connect_halfedges(left_last_hid, second_he);
        self.core_data.connect_halfedges(second_he, left_first_hid);

        // h-e
        let new_e = self.new_edges(1);
        self.he_edge_arr[first_he] = new_e;
        self.he_edge_arr[second_he] = new_e;

        // h-f
        let new_f = self.new_faces(1);
        self.he_face_arr[first_he] = fid;

        {
            let mut hid = second_he;
            loop {
                self.he_face_arr[hid] = new_f;
                hid = self.he_next(hid);
                if hid == second_he {
                    break;
                }
            }
        }

        // e-h
        self.e_halfedge_arr[new_e] = first_he;

        // f-h
        self.core_data.f_halfedge_arr[fid] = first_he;
        self.core_data.f_halfedge_arr[new_f] = second_he;
        second_he
    }

    pub fn add_face_by_halfedges<A: Allocator + Copy>(
        &mut self,
        halfedges: &[HalfedgeId],
        bump: A,
    ) -> FaceId {
        let first_new_hid = self.new_halfedges(halfedges.len());
        let new_f = self.new_faces(1);

        for (&old_hid, idx) in halfedges
            .iter()
            .zip(first_new_hid.0..self.n_halfedges_capacity())
        {
            let new_hid = idx.into();
            // h-v
            self.core_data.he_vertex_arr[new_hid] = self.core_data.he_vertex_arr[old_hid];

            // h-h
            // sibling
            insert_halfedge(&mut self.he_sibling_arr, old_hid, new_hid);
            // incoming
            insert_halfedge(&mut self.he_vert_in_next_arr, old_hid, new_hid);

            // h-e
            self.he_edge_arr[new_hid] = self.he_edge_arr[old_hid];

            // h-f
            self.he_face_arr[new_hid] = new_f;
        }

        // h-h: next halfedge
        for (ha, hb) in (first_new_hid.0..self.n_halfedges_capacity()).circular_tuple_windows() {
            self.core_data.connect_halfedges(ha.into(), hb.into());
        }

        // f-h
        self.core_data.f_halfedge_arr[new_f] = first_new_hid;
        new_f
    }

    #[inline(always)]
    pub fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.he_face_arr[hid]
    }
}

#[inline(always)]
fn insert_halfedge(next_arr: &mut [HalfedgeId], start: HalfedgeId, new_he: HalfedgeId) {
    let nnext = next_arr[start.0];
    next_arr[start] = new_he;
    next_arr[new_he] = nnext;
}

impl Mesh for SurfaceMesh {
    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        self.use_implicit_twin
    }

    #[inline(always)]
    fn core_data(&self) -> &MeshCoreData {
        &self.core_data
    }

    #[inline(always)]
    fn n_edges(&self) -> usize {
        self.n_edges
    }

    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        self.e_halfedge_arr.len()
    }

    #[inline(always)]
    fn edge_is_valid(&self, eid: EdgeId) -> bool {
        self.e_halfedge_arr[eid].valid()
    }

    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        let mut curr = hid;
        loop {
            let next = self.he_sibling_arr[curr];
            if next == hid {
                return hid;
            }

            if self.he_to(curr) != self.he_to(next) {
                return next;
            }

            curr = next;
        }
    }

    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_sibling_arr[hid]
    }

    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_vert_in_next_arr[hid]
    }

    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        self.he_edge_arr[hid]
    }

    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        self.e_halfedge_arr[eid]
    }
}
