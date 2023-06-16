use std::{cell::RefCell, collections::HashMap, rc::Weak};

use crate::{build_connect_info, mesh::Mesh, INVALID_IND};

use super::{
    ele_ranges, Edge, EdgeData, EdgeId, EdgeIter, ElementId, Face, FaceData, FaceId, FaceIter,
    FaceOrBoundaryLoopId, HalfedgeData, HalfedgeId, HalfedgeIter, MeshData, Vertex, VertexData,
    VertexId, VertexIter,
};
use bumpalo::Bump;
use itertools::Itertools;

pub struct SurfaceMesh {
    pub v_halfedge_arr: Vec<HalfedgeId>,
    pub he_next_arr: Vec<HalfedgeId>,
    pub he_vertex_arr: Vec<VertexId>,
    pub he_face_arr: Vec<FaceId>,
    pub f_halfedge_arr: Vec<HalfedgeId>,

    pub n_vertices: usize,
    pub n_halfedges: usize,
    pub n_edges: usize,
    pub n_faces: usize,

    pub he_edge_arr: Vec<EdgeId>,
    pub he_vert_in_next_arr: Vec<HalfedgeId>,
    pub he_vert_out_next_arr: Vec<HalfedgeId>,
    pub he_sibling_arr: Vec<HalfedgeId>,
    pub e_halfedge_arr: Vec<HalfedgeId>,

    pub vertices_data: Vec<Weak<RefCell<dyn MeshData<Id = VertexId> + 'static>>>,
    pub halfedges_data: Vec<Weak<RefCell<dyn MeshData<Id = HalfedgeId> + 'static>>>,
    pub edges_data: Vec<Weak<RefCell<dyn MeshData<Id = EdgeId> + 'static>>>,
    pub faces_data: Vec<Weak<RefCell<dyn MeshData<Id = FaceId> + 'static>>>,
}

impl SurfaceMesh {
    fn vertex_cycle(&self, incoming: bool) -> (Vec<HalfedgeId>, Vec<usize>) {
        let mut v_degree = vec![0usize; self.n_vertices_capacity()];
        self.halfedges().for_each(|hid| {
            let vid = if incoming {
                self.he_tip_vertex(hid)
            } else {
                self.he_vertex(hid)
            };
            v_degree[vid] += 1;
        });
        let mut vertex_separators = vec![0];
        vertex_separators.extend(v_degree.iter().scan(0, |sum, &count| {
            *sum += count;
            Some(*sum)
        }));
        let mut he_positions = vertex_separators.clone();
        let mut vertex_halfedges = vec![HalfedgeId::from(0); self.n_halfedges_capacity()];
        self.halfedges().for_each(|hid| {
            let vid = if incoming {
                self.he_tip_vertex(hid)
            } else {
                self.he_vertex(hid)
            };
            let pos = he_positions[vid];
            vertex_halfedges[pos] = hid;
            he_positions[vid] += 1;
        });
        (vertex_halfedges, vertex_separators)
    }

    #[inline]
    fn new_vertices(&mut self, n: usize) -> VertexId {
        let vid = VertexId::from(self.v_halfedge_arr.len());
        let new_len = self.n_vertices_capacity() + n;
        self.v_halfedge_arr.resize(new_len, HalfedgeId::new());
        self.n_vertices += n;
        for data in &self.vertices_data {
            if let Some(data) = data.upgrade() {
                data.borrow_mut().expand(self.n_vertices_capacity())
            }
        }
        vid
    }

    #[inline]
    fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        let hid = HalfedgeId::from(self.he_next_arr.len());
        let new_len = self.n_halfedges_capacity() + n;
        self.he_next_arr.resize(new_len, HalfedgeId::new());
        self.he_vertex_arr.resize(new_len, VertexId::new());

        self.he_face_arr.resize(new_len, FaceId::new());
        self.he_sibling_arr.resize(new_len, HalfedgeId::new());
        self.he_edge_arr.resize(new_len, EdgeId::new());
        self.he_vert_in_next_arr.resize(new_len, HalfedgeId::new());
        self.he_vert_out_next_arr.resize(new_len, HalfedgeId::new());
        self.n_halfedges += n;

        for data in &self.halfedges_data {
            if let Some(data) = data.upgrade() {
                data.borrow_mut().expand(self.n_halfedges_capacity());
            }
        }

        hid
    }

    #[inline]
    fn new_edges(&mut self, n: usize) -> EdgeId {
        let eid = self.e_halfedge_arr.len();
        let new_len = self.n_edges_capacity() + n;
        self.e_halfedge_arr.resize(new_len, HalfedgeId::new());

        for data in &self.edges_data {
            if let Some(data) = data.upgrade() {
                data.borrow_mut().expand(self.n_edges_capacity());
            }
        }
        self.n_edges += n;
        eid.into()
    }

    #[inline]
    fn new_faces(&mut self, n: usize) -> FaceId {
        let fid = self.f_halfedge_arr.len();
        let new_len = self.n_faces_capacity() + n;
        self.f_halfedge_arr.resize(new_len, HalfedgeId::new());
        for data in &self.faces_data {
            if let Some(data) = data.upgrade() {
                data.borrow_mut().expand(self.n_faces_capacity());
            }
        }
        self.n_faces += n;
        fid.into()
    }

    #[inline]
    pub fn split_edge(&mut self, eid: EdgeId, bump: &Bump) -> VertexId {
        let mut e_halfedges = Vec::new_in(bump);
        e_halfedges.extend(self.edge(eid).halfedges());
        let (va, vb) = (
            self.he_vertex(e_halfedges[0]),
            self.he_tip_vertex(e_halfedges[0]),
        );
        let new_v = self.new_vertices(1);
        let new_halfedges = ele_ranges::<HalfedgeId>(
            self.new_halfedges(e_halfedges.len()).0,
            e_halfedges.len(),
            bump,
        );

        // v-h
        self.v_halfedge_arr[new_v] = new_halfedges[0];

        // h-v
        for &hid in &new_halfedges {
            self.he_vertex_arr[hid] = new_v;
        }

        // h-h
        let he_pos_map = HashMap::<HalfedgeId, usize>::from_iter(
            e_halfedges.iter().enumerate().map(|(i, &hid)| (hid, i)),
        );
        let to_new_halfedge = |hid: HalfedgeId, positions: &mut Vec<usize, _>| -> HalfedgeId {
            if let Some(&pos) = he_pos_map.get(&hid) {
                positions.push(pos);
                new_halfedges[pos]
            } else {
                hid
            }
        };
        let mut e_va_in_positions = Vec::new_in(bump);
        let mut e_vb_in_positions = Vec::new_in(bump);
        let mut va_in_halfedges = Vec::new_in(bump);
        va_in_halfedges.extend(
            self.vertex(va)
                .incoming_halfedge()
                .map(|hid| to_new_halfedge(hid, &mut e_va_in_positions)),
        );
        let mut vb_in_halfedges = Vec::new_in(bump);
        vb_in_halfedges.extend(
            self.vertex(vb)
                .incoming_halfedge()
                .map(|hid| to_new_halfedge(hid, &mut e_vb_in_positions)),
        );
        for (&ha, &hb) in va_in_halfedges.iter().circular_tuple_windows() {
            self.he_vert_in_next_arr[ha] = hb;
        }
        for (&ha, &hb) in vb_in_halfedges.iter().circular_tuple_windows() {
            self.he_vert_in_next_arr[ha] = hb;
        }
        for (&ha, &hb) in e_halfedges.iter().circular_tuple_windows() {
            self.he_vert_in_next_arr[ha] = hb;
        }
        for (&ha, &hb) in new_halfedges.iter().circular_tuple_windows() {
            self.he_vert_out_next_arr[ha] = hb;
            // self.he_sibling_arr[ha] = hb;
        }

        let mut e_va_halfedges = Vec::new_in(bump);
        e_va_halfedges.extend(
            e_vb_in_positions
                .iter()
                .map(|&idx| e_halfedges[idx])
                .chain(e_va_in_positions.iter().map(|&idx| new_halfedges[idx])),
        );
        let mut e_vb_halfedges = Vec::new_in(bump);
        e_vb_halfedges.extend(
            e_va_in_positions
                .into_iter()
                .map(|idx| e_halfedges[idx])
                .chain(e_vb_in_positions.into_iter().map(|idx| new_halfedges[idx])),
        );

        for (&ha, &hb) in e_va_halfedges.iter().circular_tuple_windows() {
            self.he_sibling_arr[ha] = hb;
        }
        for (&ha, &hb) in e_vb_halfedges.iter().circular_tuple_windows() {
            self.he_sibling_arr[ha] = hb;
        }

        for (&prev, &next) in e_halfedges.iter().zip(&new_halfedges) {
            insert_halfedge(&mut self.he_next_arr, prev, next);
        }

        // h-e
        for &he in &e_va_halfedges {
            self.he_edge_arr[he] = eid;
        }
        let new_e = self.new_edges(1);
        for &he in &e_vb_halfedges {
            self.he_edge_arr[he] = new_e;
        }

        // e-h
        self.e_halfedge_arr[eid] = e_va_halfedges[0];
        self.e_halfedge_arr[new_e] = e_vb_halfedges[0];

        // h-f
        for (old_he, new_he) in e_halfedges.into_iter().zip(new_halfedges) {
            self.he_face_arr[new_he] = self.he_face_arr[old_he];
        }

        new_v
    }

    #[inline]
    pub fn split_face(
        &mut self,
        fid: FaceId,
        va: VertexId,
        vb: VertexId,
        bump: &Bump,
    ) -> HalfedgeId {
        let mut face_halfedges = Vec::new_in(bump);
        face_halfedges.extend(self.face(fid).halfedges());
        let (mut pos1, mut pos2) = (INVALID_IND, INVALID_IND);
        for (i, &hid) in face_halfedges.iter().enumerate() {
            let v = self.he_vertex(hid);
            if v == va || v == vb {
                if pos1 == INVALID_IND {
                    pos1 = i;
                } else {
                    pos2 = i;
                    break;
                }
            }
        }
        let (va, vb) = if self.he_vertex(face_halfedges[pos1]) == va {
            (va, vb)
        } else {
            (vb, va)
        };
        let mut first_halfedges = Vec::new_in(bump);
        first_halfedges.extend(face_halfedges[pos1..pos2].iter().map(|&hid| hid));
        let mut second_halfedges = Vec::new_in(bump);
        second_halfedges.extend(
            face_halfedges[pos2..]
                .iter()
                .chain(&face_halfedges[..pos1])
                .map(|&hid| hid),
        );
        let first_he = self.new_halfedges(2);
        let second_he = HalfedgeId::from(first_he.0 + 1);

        // h-v
        self.he_vertex_arr[first_he] = vb;
        self.he_vertex_arr[second_he] = va;

        // h-h
        // sibling halfedges
        self.he_sibling_arr[first_he] = second_he;
        self.he_sibling_arr[second_he] = first_he;
        // incoming halfedges and outcoming halfedges
        insert_halfedge(
            &mut self.he_vert_in_next_arr,
            *second_halfedges.last().unwrap(),
            first_he,
        );
        insert_halfedge(
            &mut self.he_vert_out_next_arr,
            second_halfedges[0],
            first_he,
        );
        insert_halfedge(
            &mut self.he_vert_in_next_arr,
            *first_halfedges.last().unwrap(),
            second_he,
        );
        insert_halfedge(
            &mut self.he_vert_out_next_arr,
            first_halfedges[0],
            second_he,
        );
        // next halfedges
        self.he_next_arr[first_he] = first_halfedges[0];
        self.he_next_arr[*first_halfedges.last().unwrap()] = first_he;
        self.he_next_arr[second_he] = second_halfedges[0];
        self.he_next_arr[*second_halfedges.last().unwrap()] = second_he;

        // h-e
        let new_e = self.new_edges(1);
        self.he_edge_arr[first_he] = new_e;
        self.he_edge_arr[second_he] = new_e;

        // h-f
        let new_f = self.new_faces(1);
        self.he_face_arr[first_he] = fid;
        self.he_face_arr[second_he] = new_f;
        for hid in second_halfedges {
            self.he_face_arr[hid] = new_f;
        }

        // e-h
        self.e_halfedge_arr[new_e] = first_he;

        // f-h
        self.f_halfedge_arr[new_f] = second_he;
        second_he
    }
}

#[inline(always)]
fn insert_halfedge(next_arr: &mut [HalfedgeId], start: HalfedgeId, new_he: HalfedgeId) {
    let nnext = next_arr[start.0];
    next_arr[start.0] = new_he;
    next_arr[new_he.0] = nnext;
}

impl From<Vec<Vec<usize>>> for SurfaceMesh {
    fn from(polygons: Vec<Vec<usize>>) -> Self {
        let n_faces = polygons.len();
        let mut n_vertices = 0usize;
        polygons.iter().for_each(|polygon| {
            polygon
                .iter()
                .for_each(|vid| n_vertices = n_vertices.max(*vid));
        });
        n_vertices += 1;
        let mut mesh = Self {
            v_halfedge_arr: vec![HalfedgeId::new(); n_vertices],
            he_next_arr: Vec::new(),
            he_vertex_arr: Vec::new(),
            he_face_arr: Vec::new(),
            f_halfedge_arr: vec![HalfedgeId::new(); n_faces],

            n_vertices,
            n_halfedges: 0,
            n_edges: 0,
            n_faces,

            he_edge_arr: Vec::new(),
            he_vert_in_next_arr: Vec::new(),
            he_vert_out_next_arr: Vec::new(),
            he_sibling_arr: Vec::new(),
            e_halfedge_arr: Vec::new(),

            vertices_data: Vec::new(),
            halfedges_data: Vec::new(),
            edges_data: Vec::new(),
            faces_data: Vec::new(),
        };

        for (fid, polygon) in polygons.iter().enumerate() {
            let mut first_hid = INVALID_IND.into();
            let mut prev_hid = first_hid;
            for i in 0..polygon.len() {
                let tail = polygon[i];
                let hid = mesh.new_halfedges(1);

                mesh.he_vertex_arr[hid] = tail.into();
                mesh.he_face_arr[hid] = fid.into();

                mesh.v_halfedge_arr[tail] = hid;
                if i == 0 {
                    mesh.f_halfedge_arr[fid] = hid;
                    first_hid = hid;
                } else {
                    mesh.he_next_arr[prev_hid] = hid;
                }
                prev_hid = hid;
            }
            mesh.he_next_arr[prev_hid] = first_hid;
        }

        let mut edge_history = HashMap::<(usize, usize), usize>::new();

        // build edge
        {
            let mut hid = 0;
            for polygon in polygons {
                for i in 0..polygon.len() {
                    let tail = polygon[i];
                    let tip = polygon[(i + 1) % polygon.len()];
                    let key = if tail < tip { (tail, tip) } else { (tip, tail) };
                    if let Some(prev_hid) = edge_history.get_mut(&key) {
                        // We're already seen this edge, connect to the previous halfedge incident on the edge
                        mesh.he_sibling_arr[hid] = (*prev_hid).into();
                        let eid = mesh.he_edge_arr[*prev_hid];
                        mesh.he_edge_arr[hid] = eid;
                        *prev_hid = hid;
                    } else {
                        // This is the first time we've ever seen this edge, create a new edge object
                        let new_eid = mesh.new_edges(1);
                        mesh.he_edge_arr[hid] = new_eid;
                        mesh.he_sibling_arr[hid] = HalfedgeId::new();
                        mesh.e_halfedge_arr[new_eid] = hid.into();
                        edge_history.insert(key, hid);
                    }
                    hid += 1;
                }
            }
        }
        // Complete the sibling cycle by following backwards each edge until we reach the first sibling-less entry
        for last_he in edge_history.into_values() {
            if !mesh.he_sibling_arr[last_he].valid() {
                // Any edges which never got any sibling entries at all are boundary halfedges
                mesh.he_sibling_arr[last_he] = last_he.into();
                continue;
            }

            // Get the index of the first halfedge in the sibling cycle to complete the cycle
            let mut curr_he = last_he;
            while mesh.he_sibling_arr[curr_he].valid() {
                curr_he = *mesh.he_sibling_arr[curr_he];
            }
            mesh.he_sibling_arr[curr_he] = last_he.into(); // connect the first to the last
        }

        let (v_in_halfedges, v_in_separators) = mesh.vertex_cycle(true);
        let (v_out_halfedges, v_out_separators) = mesh.vertex_cycle(false);
        let n_vertices = mesh.n_vertices();
        for idx in 0..n_vertices {
            let vid = VertexId::from(idx);
            if !mesh.vertex_is_valid(vid) {
                continue;
            }
            {
                let (start, end) = (v_in_separators[vid], v_in_separators[*vid + 1]);
                let len = end - start;
                for i in start..end {
                    let ha = v_in_halfedges[i];
                    let hb = v_in_halfedges[start + (i - start + 1) % len];
                    mesh.he_vert_in_next_arr[ha] = hb;
                }
            }
            {
                let (start, end) = (v_out_separators[vid], v_out_separators[*vid + 1]);
                let len = end - start;
                for i in start..end {
                    let ha = v_out_halfedges[i];
                    let hb = v_out_halfedges[start + (i - start + 1) % len];
                    mesh.he_vert_out_next_arr[ha] = hb;
                }
            }
        }
        mesh
    }
}

impl Mesh for SurfaceMesh {
    build_connect_info!();

    /// the number of edges
    #[inline(always)]
    fn n_edges(&self) -> usize {
        return self.n_edges;
    }

    /// the capacity of edges
    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        return self.e_halfedge_arr.len();
    }

    /// edge id is valid
    #[inline(always)]
    fn edge_is_valid(&self, eid: super::EdgeId) -> bool {
        self.e_halfedge_arr[eid].valid()
    }

    /// the twin halfedge of the halfedge
    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_sibling_arr[hid])
    }

    /// the sibling halfedge of the halfedge
    #[inline(always)]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_sibling_arr[hid])
    }

    /// the next incoming halfedge of the halfedge
    #[inline(always)]
    fn he_next_incoming_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_vert_in_next_arr[hid]
    }

    /// the next outgoing halfedge of the halfedge
    #[inline(always)]
    fn he_next_outgoing_neighbor(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_vert_out_next_arr[hid])
    }

    /// the edge of this halfedge
    #[inline(always)]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        EdgeId::from(self.he_edge_arr[hid])
    }

    /// the face or boundary loop of the halfedge
    #[inline(always)]
    fn he_face_or_boundary_loop(&self, hid: HalfedgeId) -> FaceOrBoundaryLoopId {
        FaceOrBoundaryLoopId::Face(self.he_face_arr[hid])
    }

    /// the first halfedge of the edge
    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        HalfedgeId::from(self.e_halfedge_arr[eid])
    }

    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        false
    }
}
