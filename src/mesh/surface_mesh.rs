use std::{cell::RefCell, collections::HashMap, rc::Weak};

use crate::{build_connect_info, mesh::Mesh, INVALID_IND};

use super::{
    EdgeId, EdgeIter, ElementId, FaceId, FaceIter, FaceOrBoundaryLoopId, HalfedgeData, HalfedgeId,
    HalfedgeIter, MeshData, VertexId, VertexIter,
};
use bumpalo::{collections::Vec, Bump};

pub struct SurfaceMesh<'b> {
    v_halfedge_arr: Vec<'b, HalfedgeId>,
    he_next_arr: Vec<'b, HalfedgeId>,
    he_vertex_arr: Vec<'b, VertexId>,
    he_face_arr: Vec<'b, FaceId>,
    f_halfedge_arr: Vec<'b, HalfedgeId>,

    n_vertices: usize,
    n_halfedges: usize,
    n_edges: usize,
    n_faces: usize,

    he_edge_arr: Vec<'b, EdgeId>,
    he_vert_in_next_arr: Vec<'b, HalfedgeId>,
    he_vert_in_prev_arr: Vec<'b, HalfedgeId>,
    he_vert_out_next_arr: Vec<'b, HalfedgeId>,
    he_vert_out_prev_arr: Vec<'b, HalfedgeId>,
    he_sibling_arr: Vec<'b, HalfedgeId>,
    e_halfedge_arr: Vec<'b, HalfedgeId>,

    pub halfedges_data: std::vec::Vec<Weak<RefCell<dyn MeshData<Id = HalfedgeId> + 'b>>>,
}

impl<'b> SurfaceMesh<'b> {
    fn vertex_cycle(&self, incoming: bool) -> (Vec<'b, HalfedgeId>, Vec<'b, usize>) {
        let bump = self.v_halfedge_arr.bump();
        let mut v_degree = bumpalo::vec![in &bump; 0usize; self.n_vertices_capacity()];
        self.halfedges().for_each(|hid| {
            let vid = if incoming {
                self.he_tip_vertex(hid)
            } else {
                self.he_vertex(hid)
            };
            v_degree[vid] += 1;
        });
        let mut vertex_separators = bumpalo::vec![in bump; 0];
        vertex_separators.extend(v_degree.iter().scan(0, |sum, &count| {
            *sum += count;
            Some(*sum)
        }));
        let mut he_positions = vertex_separators.clone();
        let mut vertex_halfedges =
            bumpalo::vec![in bump; HalfedgeId::from(0); self.n_halfedges_capacity()];
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

    fn new_halfedge(&mut self) -> HalfedgeId {
        let hid = HalfedgeId::from(self.he_next_arr.len());
        self.he_next_arr.push(HalfedgeId::new());
        self.he_vertex_arr.push(VertexId::new());

        self.he_face_arr.push(FaceId::new());
        self.he_sibling_arr.push(HalfedgeId::new());
        self.he_edge_arr.push(EdgeId::new());
        self.he_vert_in_next_arr.push(HalfedgeId::new());
        self.he_vert_in_prev_arr.push(HalfedgeId::new());
        self.he_vert_out_next_arr.push(HalfedgeId::new());
        self.he_vert_out_prev_arr.push(HalfedgeId::new());
        self.n_halfedges += 1;

        for data in &self.halfedges_data {
            if let Some(data) = data.upgrade() {
                data.borrow_mut().expand(self.n_halfedges);
            }
        }

        hid
    }

    fn new_edge(&mut self) -> EdgeId {
        let eid = self.e_halfedge_arr.len();
        self.e_halfedge_arr.push(HalfedgeId::new());

        // TODO: update edge callback

        self.n_edges += 1;
        eid.into()
    }
}

impl<'b> From<Vec<'b, Vec<'b, usize>>> for SurfaceMesh<'b> {
    fn from(polygons: Vec<'b, Vec<'b, usize>>) -> Self {
        let n_faces = polygons.len();
        let mut n_vertices = 0usize;
        polygons.iter().for_each(|polygon| {
            polygon
                .iter()
                .for_each(|vid| n_vertices = n_vertices.max(*vid));
        });
        let bump = polygons.bump();
        n_vertices += 1;
        let mut mesh = Self {
            v_halfedge_arr: bumpalo::vec![in bump; HalfedgeId::new(); n_vertices],
            he_next_arr: Vec::new_in(bump),
            he_vertex_arr: Vec::new_in(bump),
            he_face_arr: Vec::new_in(bump),
            f_halfedge_arr: bumpalo::vec![in bump; HalfedgeId::new(); n_faces],

            n_vertices,
            n_halfedges: 0,
            n_edges: 0,
            n_faces,

            he_edge_arr: Vec::new_in(bump),
            he_vert_in_next_arr: Vec::new_in(bump),
            he_vert_in_prev_arr: Vec::new_in(bump),
            he_vert_out_next_arr: Vec::new_in(bump),
            he_vert_out_prev_arr: Vec::new_in(bump),
            he_sibling_arr: Vec::new_in(bump),
            e_halfedge_arr: Vec::new_in(bump),

            halfedges_data: std::vec::Vec::new(),
        };

        for (fid, polygon) in polygons.iter().enumerate() {
            let mut first_hid = INVALID_IND.into();
            let mut prev_hid = first_hid;
            for i in 0..polygon.len() {
                let tail = polygon[i];
                let hid = mesh.new_halfedge();

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
                        let new_eid = mesh.new_edge();
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
                    mesh.he_vert_in_prev_arr[hb] = ha;
                }
            }
            {
                let (start, end) = (v_out_separators[vid], v_out_separators[*vid + 1]);
                let len = end - start;
                for i in start..end {
                    let ha = v_out_halfedges[i];
                    let hb = v_out_halfedges[start + (i - start + 1) % len];
                    mesh.he_vert_out_next_arr[ha] = hb;
                    mesh.he_vert_out_prev_arr[hb] = ha;
                }
            }
        }
        mesh
    }
}

impl<'b> Mesh<'b> for SurfaceMesh<'b> {
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
