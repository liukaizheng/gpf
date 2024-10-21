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
                mesh.core_data.he_vertex_arr[hid] = b.into();
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
    #[inline(always)]
    pub fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.he_face_arr[hid]
    }
    #[inline(always)]
    pub fn he_is_boundary(&self, hid: HalfedgeId) -> bool {
        !self.he_face(hid).valid()
    }

    #[inline(always)]
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
    fn edge_is_valid(&self, eid: EdgeId) -> bool {
        let idx = eid.0 << 1;
        self.halfedge_is_valid(idx.into()) || self.halfedge_is_valid((idx + 1).into())
    }

    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
        let mut curr = hid;
        loop {
            let next = self.he_next(curr);
            if next == hid {
                return curr;
            }
            curr = next;
        }
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
        self.he_twin(self.he_next(hid))
    }

    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        (hid.0 >> 1).into()
    }

    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        let idx = eid.0 << 1;
        if self.halfedge_is_valid(idx.into()) {
            idx.into()
        } else {
            (idx + 1).into()
        }
    }
}
