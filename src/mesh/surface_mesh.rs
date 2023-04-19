use std::collections::HashMap;

use crate::{build_connect_info, mesh::Mesh, INVALID_IND};

use super::{
    BoundaryLoopId, BoundaryLoopIter, EdgeId, EdgeIter, FaceId, FaceIter, FaceOrBoundaryLoopId,
    HalfedgeId, HalfedgeIter, VertexId, VertexIter, BL_START,
};

struct SurfaceMesh {
    v_halfedge_arr: Vec<usize>,
    he_next_arr: Vec<usize>,
    he_vertex_arr: Vec<usize>,
    he_face_arr: Vec<usize>,
    f_halfedge_arr: Vec<usize>,
    bl_halfedge_arr: Vec<usize>,

    n_vertices: usize,
    n_halfedges: usize,
    n_interior_halfedges: usize,
    n_edges: usize,
    n_faces: usize,
    n_loops: usize,

    v_he_in_start_arr: Vec<usize>,
    v_he_out_start_arr: Vec<usize>,
    he_edge_arr: Vec<usize>,
    he_vert_in_next_arr: Vec<usize>,
    he_vert_in_prev_arr: Vec<usize>,
    he_vert_out_next_arr: Vec<usize>,
    he_vert_out_prev_arr: Vec<usize>,
    he_sibling_arr: Vec<usize>,
    e_halfedge_arr: Vec<usize>,
}

impl SurfaceMesh {
    fn generate_vertex_iter(&self, incoming: bool) -> (Vec<HalfedgeId>, Vec<usize>) {
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
            v_halfedge_arr: vec![INVALID_IND; n_vertices],
            he_next_arr: vec![],
            he_vertex_arr: vec![],
            he_face_arr: vec![],
            f_halfedge_arr: vec![INVALID_IND; n_faces],
            bl_halfedge_arr: vec![],

            n_vertices,
            n_halfedges: 0,
            n_interior_halfedges: 0,
            n_edges: 0,
            n_faces,
            n_loops: 0,

            v_he_in_start_arr: vec![],
            v_he_out_start_arr: vec![],
            he_edge_arr: vec![],
            he_vert_in_next_arr: vec![],
            he_vert_in_prev_arr: vec![],
            he_vert_out_next_arr: vec![],
            he_vert_out_prev_arr: vec![],
            he_sibling_arr: vec![],
            e_halfedge_arr: vec![],
        };

        for (fid, polygon) in polygons.iter().enumerate() {
            let mut first_hid = INVALID_IND.into();
            let mut prev_hid = first_hid;
            for i in 0..polygon.len() {
                let tail = polygon[i];
                let hid = mesh.new_halfedge(true);

                mesh.he_vertex_arr[hid] = tail;
                mesh.he_face_arr[hid] = fid;

                mesh.v_halfedge_arr[tail] = hid.0;
                if i == 0 {
                    mesh.f_halfedge_arr[fid] = hid.0;
                    first_hid = hid;
                } else {
                    mesh.he_next_arr[prev_hid] = hid.0;
                }
                prev_hid = hid;
            }
            mesh.he_next_arr[prev_hid] = first_hid.0;
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
                        mesh.he_sibling_arr[hid] = *prev_hid;
                        let eid = mesh.he_edge_arr[*prev_hid];
                        mesh.he_edge_arr[hid] = eid;
                        *prev_hid = hid;
                    } else {
                        // This is the first time we've ever seen this edge, create a new edge object
                        let new_eid = mesh.new_edge();
                        mesh.he_edge_arr[hid] = new_eid.0;
                        mesh.he_sibling_arr[hid] = INVALID_IND;
                        mesh.e_halfedge_arr[new_eid] = hid;
                        edge_history.insert(key, hid);
                    }
                    hid += 1;
                }
            }
        }
        // Complete the sibling cycle by follwing backwards each edge until we reach the first sibling-less entry
        for last_he in edge_history.into_values() {
            if mesh.he_sibling_arr[last_he] == INVALID_IND {
                // Any edges which never got any sibling entries at all are boundary halfedges
                mesh.he_sibling_arr[last_he] = last_he;
                continue;
            }

            // Get the index of the first halfedge in the sibling cycle to complete the cycle
            let mut curr_he = last_he;
            while mesh.he_sibling_arr[curr_he] != INVALID_IND {
                curr_he = mesh.he_sibling_arr[curr_he];
            }
            mesh.he_sibling_arr[curr_he] = last_he; // connect the first to the last
        }
        mesh
    }
}

impl Mesh for SurfaceMesh {
    build_connect_info!();

    /// the number of halfedges
    #[inline(always)]
    fn n_halfedges(&self) -> usize {
        return self.n_halfedges;
    }

    /// the number of edges
    #[inline(always)]
    fn n_edges(&self) -> usize {
        return self.n_edges;
    }

    /// the number of faces
    #[inline(always)]
    fn n_faces(&self) -> usize {
        return self.n_faces;
    }

    /// the number of boundary loops
    #[inline(always)]
    fn n_boundary_loops(&self) -> usize {
        return self.n_loops;
    }

    /// the capacity of vertices
    #[inline(always)]
    fn n_vertices_capacity(&self) -> usize {
        return self.v_halfedge_arr.len();
    }

    /// the capacity of halfedges
    #[inline(always)]
    fn n_halfedges_capacity(&self) -> usize {
        return self.he_next_arr.len();
    }

    /// the capacity of edges
    #[inline(always)]
    fn n_edges_capacity(&self) -> usize {
        return self.e_halfedge_arr.len();
    }

    /// the capacity of faces
    #[inline(always)]
    fn n_faces_capacity(&self) -> usize {
        return self.f_halfedge_arr.len();
    }

    /// the capacity of boundary loops
    #[inline(always)]
    fn n_boundary_loop_capacity(&self) -> usize {
        return self.bl_halfedge_arr.len();
    }

    /// vertex id is valid
    #[inline(always)]
    fn vertex_is_valid(&self, vid: VertexId) -> bool {
        return self.v_halfedge_arr[vid] != INVALID_IND;
    }

    /// halfedge id is valid
    #[inline(always)]
    fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool {
        return self.he_next_arr[hid] != INVALID_IND;
    }

    /// edge id is valid
    #[inline(always)]
    fn edge_is_valid(&self, eid: super::EdgeId) -> bool {
        return self.e_halfedge_arr[eid] != INVALID_IND;
    }

    /// face id is valid
    #[inline(always)]
    fn face_is_valid(&self, fid: FaceId) -> bool {
        return self.f_halfedge_arr[fid] != INVALID_IND;
    }

    /// boundary loop id is valid
    #[inline(always)]
    fn boundary_loop_is_valid(&self, blid: BoundaryLoopId) -> bool {
        return self.bl_halfedge_arr[blid] != INVALID_IND;
    }

    /// start vertex iterator
    #[inline(always)]
    fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self> {
        VertexIter::new(vid, self)
    }

    /// start halfedge iterator
    #[inline(always)]
    fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self> {
        HalfedgeIter::new(hid, self)
    }

    /// start edge iterator
    #[inline(always)]
    fn edge<'a>(&'a self, eid: EdgeId) -> super::EdgeIter<'a, Self> {
        EdgeIter::new(eid, self)
    }

    /// start face iterator
    #[inline(always)]
    fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self> {
        FaceIter::new(fid, self)
    }

    /// start boundary loop iterator
    #[inline(always)]
    fn boundary_loop<'a>(&'a self, blid: BoundaryLoopId) -> BoundaryLoopIter<'a, Self> {
        BoundaryLoopIter::new(blid, self)
    }

    #[inline(always)]
    fn vertices(&self) -> VertexIter<'_, Self> {
        VertexIter::new(VertexId::from(0), self)
    }

    #[inline(always)]
    fn halfedges(&self) -> HalfedgeIter<'_, Self> {
        HalfedgeIter::new(HalfedgeId::from(0), self)
    }

    #[inline(always)]
    fn edges(&self) -> EdgeIter<'_, Self> {
        EdgeIter::new(EdgeId::from(0), self)
    }

    #[inline(always)]
    fn faces(&self) -> FaceIter<'_, Self> {
        FaceIter::new(FaceId::from(0), self)
    }

    /// the halfedge starting from this vertex
    #[inline(always)]
    fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
        HalfedgeId::from(self.v_halfedge_arr[vid])
    }

    /// the start vertex of the halfedge
    #[inline(always)]
    fn he_vertex(&self, hid: HalfedgeId) -> VertexId {
        VertexId::from(self.he_vertex_arr[hid])
    }

    /// the end vertex of the halfedge
    #[inline(always)]
    fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId {
        VertexId::from(self.he_vertex_arr[self.he_next_arr[hid]])
    }

    /// the next halfedge of the halfedge
    #[inline(always)]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        HalfedgeId::from(self.he_next_arr[hid])
    }

    /// the previous halfedge of the halfedge
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
        HalfedgeId::from(self.he_vert_in_next_arr[hid])
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
        let id = self.he_face_arr[hid];
        if id == INVALID_IND {
            FaceOrBoundaryLoopId::INVALID
        } else if id < self.he_face_arr.len() {
            FaceOrBoundaryLoopId::Face(FaceId::from(id))
        } else {
            FaceOrBoundaryLoopId::BoundaryLoop(BoundaryLoopId::from(BL_START - id))
        }
    }

    /// the first halfedge of the edge
    #[inline(always)]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        HalfedgeId::from(self.e_halfedge_arr[eid])
    }

    /// the first halfedge of the face
    #[inline(always)]
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
        HalfedgeId::from(self.f_halfedge_arr[fid])
    }

    #[inline(always)]
    fn use_implicit_twin(&self) -> bool {
        true
    }

    fn new_halfedge(&mut self, is_interior: bool) -> HalfedgeId {
        let hid = HalfedgeId::from(self.he_next_arr.len());
        self.he_next_arr.push(INVALID_IND);
        self.he_vertex_arr.push(INVALID_IND);

        self.he_face_arr.push(INVALID_IND);
        self.he_sibling_arr.push(INVALID_IND);
        self.he_edge_arr.push(INVALID_IND);
        self.he_vert_in_next_arr.push(INVALID_IND);
        self.he_vert_in_prev_arr.push(INVALID_IND);
        self.he_vert_out_next_arr.push(INVALID_IND);
        self.he_vert_out_prev_arr.push(INVALID_IND);
        self.n_halfedges += 1;

        if is_interior {
            self.n_interior_halfedges += 1;
        }

        // TODO: update callback

        hid
    }

    fn new_edge(&mut self) -> EdgeId {
        let eid = self.e_halfedge_arr.len();
        self.e_halfedge_arr.push(INVALID_IND);
        // TODO: update edge callback
        self.n_edges += 1;
        eid.into()
    }
}
