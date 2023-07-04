use std::{
    alloc::Allocator,
    time::{Duration, Instant},
};

use bumpalo::Bump;
use hashbrown::HashMap;
use itertools::Itertools;
use serde::{Serialize, Deserialize};

use crate::{
    graphcut::GraphCut,
    math::{cross, norm, sub},
    mesh::{EdgeId, Face, FaceId, HalfedgeId, Mesh, SurfaceMesh, VertexId},
    predicates::{
        orient3d::orient3d, sign_reversed, ExplicitPoint3D, ImplicitPoint3D, ImplicitPointLPI,
        ImplicitPointTPI, Orientation, Point3D,
    },
    triangle::TetMesh,
    INVALID_IND,
};

use super::{conforming_mesh::Constraints, point};

#[derive(Clone)]
struct BSPEdgeData {
    parents: Vec<VertexId>,
}

impl BSPEdgeData {
    fn new(parents: Vec<VertexId>) -> Self {
        Self { parents }
    }
}

#[derive(Serialize, Deserialize)]
pub struct Edge {
    pub e: [usize; 2],
    pub w: f64,
}

impl Edge {
    pub fn new(a: usize, b: usize, w: f64) -> Self {
        Self { e: [a, b], w }
    }
}

#[derive(Serialize, Deserialize)]
pub struct G {
    pub external: Vec<f64>,
    pub internal: Vec<f64>,
    pub edges: Vec<Edge>,
}

impl G {
    pub fn new() -> Self {
        Self {
            external: Vec::new(),
            internal: Vec::new(),
            edges: Vec::new(),
        }
    }
}

#[derive(Clone)]
struct BSPFaceData {
    cells: [usize; 2],
    // coplanar triangles
    triangles: Vec<usize>,
    plane: [VertexId; 3],
}

impl BSPFaceData {
    fn new(c1: usize, c2: usize, triangles: Vec<usize>, plane: [VertexId; 3]) -> Self {
        Self {
            cells: [c1, c2],
            triangles,
            plane,
        }
    }
}

#[derive(Clone)]
struct BSPCellData {
    faces: Vec<FaceId>,
    inner_triangles: Vec<usize>,
}

impl BSPCellData {
    fn new(faces: Vec<FaceId>, inner_triangles: Vec<usize>) -> Self {
        Self {
            faces,
            inner_triangles,
        }
    }
}

pub(crate) struct BSPComplex {
    points: Vec<Point3D>,
    pub mesh: SurfaceMesh,
    edge_data: Vec<BSPEdgeData>,
    face_data: Vec<BSPFaceData>,
    cell_data: Vec<BSPCellData>,
    constraints: Vec<VertexId>,
    n_ori_triangles: usize,

    vert_orientations: Vec<HashMap<VertexId, Orientation>>,
    vert_visits: Vec<bool>,

    edge_visits: Vec<bool>,

    pub ori_duration: Duration,
    pub split_duration: Duration,
}

#[inline(always)]
fn tet_face_is_new(tid: usize, adj_tid: usize, adj_cid: usize) -> bool {
    adj_tid > tid || adj_cid == INVALID_IND
}

impl BSPComplex {
    pub(crate) fn new(
        tet_mesh: TetMesh,
        constraints_data: &Constraints,
        tet_mark: [Vec<Vec<usize>>; 5],
    ) -> Self {
        let points = Vec::from_iter(
            tet_mesh
                .points
                .chunks(3)
                .map(|data| Point3D::Explicit(ExplicitPoint3D::from(data))),
        );

        let new_tet_orders = remove_ghost_tets(&tet_mesh);
        let mut face_positions = vec![INVALID_IND; tet_mesh.tets.len() << 2];
        let mut faces: Vec<Vec<usize>> = Vec::new();
        let mut face_data: Vec<BSPFaceData> = Vec::new();
        let mut cell_data: Vec<BSPCellData> = Vec::new();
        for (tid, tet) in tet_mesh.tets.iter().enumerate() {
            let cell_idx = new_tet_orders[tid];
            if cell_idx == INVALID_IND {
                continue;
            }

            let mut cell_faces = Vec::new();
            for (i, nei) in tet.nei.iter().enumerate() {
                let adj_cell_idx = new_tet_orders[nei.tet];
                if tet_face_is_new(tid, nei.tet, adj_cell_idx) {
                    let tri = vec![tet_mesh.org(nei), tet_mesh.apex(nei), tet_mesh.dest(nei)];
                    let mut face = BSPFaceData::new(
                        cell_idx,
                        adj_cell_idx,
                        Vec::new(),
                        [tri[0].into(), tri[1].into(), tri[2].into()],
                    );
                    faces.push(tri);
                    insert_coplanar_triangles(
                        &tet_mark[i][tid],
                        &mut face.triangles,
                        constraints_data.n_ori_triangles,
                    );
                    let fid = face_data.len();
                    face_positions[(tid << 2) + i] = fid;
                    cell_faces.push(fid.into());
                    face_data.push(face);
                } else {
                    let j = nei.ver & 3;
                    let fid = face_positions[(nei.tet << 2) + j];
                    let face = &mut face_data[fid];
                    insert_coplanar_triangles(
                        &tet_mark[j][nei.tet],
                        &mut face.triangles,
                        constraints_data.n_ori_triangles,
                    );
                    cell_faces.push(fid.into());
                }
            }
            cell_data.push(BSPCellData::new(cell_faces, tet_mark[4][tid].clone()));
        }

        let mesh = SurfaceMesh::from(faces);

        let vert_orientations = vec![HashMap::new(); constraints_data.triangles.len() / 3];
        let vert_visits = vec![false; points.len()];

        let edge_visits = vec![false; mesh.n_edges()];
        let edge_data = Vec::from_iter(
            mesh.edges()
                .map(|eid| BSPEdgeData::new(Vec::from_iter(mesh.e_vertices(eid)))),
        );

        let constraints = Vec::from_iter(constraints_data.triangles.iter().map(|&idx| idx.into()));

        Self {
            points,
            mesh,
            edge_data,
            face_data,
            cell_data,
            constraints,
            n_ori_triangles: constraints_data.n_ori_triangles,
            edge_visits,
            vert_orientations,
            vert_visits,

            ori_duration: Duration::from_millis(0),
            split_duration: Duration::from_millis(0),
        }
    }

    #[inline(always)]
    pub(crate) fn n_cells(&self) -> usize {
        self.cell_data.len()
    }

    #[inline(always)]
    pub(crate) fn splittable(&self, cid: usize) -> bool {
        !self.cell_data[cid].inner_triangles.is_empty()
    }

    pub(crate) fn split_cell<A: Allocator + Copy>(&mut self, cid: usize, bump: A) {
        let tid = *self.cell_data[cid].inner_triangles.last().unwrap();
        let tri = self.triangle(tid).to_vec_in(bump);
        let start = Instant::now();
        self.cell_data[cid].inner_triangles.pop();
        let mut coplanar_triangles = self.separate_out_coplanar_triangles(tid, cid, bump);
        if !self.is_virtual(tid) {
            coplanar_triangles.push(tid);
        }
        let (cell_verts, cell_edges) = self.cell_verts_and_edges(cid, bump);
        verts_orient_wrt_plane(
            &self.points[tri[0]],
            &self.points[tri[1]],
            &self.points[tri[2]],
            &cell_verts,
            &self.points,
            &mut self.vert_orientations[tid],
            bump,
        );
        self.ori_duration += start.elapsed();
        let start = Instant::now();

        let (mut n_over, mut n_under) = (0, 0);
        let mut cell_orientations = Vec::new();
        for vid in &cell_verts {
            cell_orientations.push(*self.vert_orientations[tid].get(vid).unwrap());
            match self.vert_orientations[tid].get(vid).unwrap() {
                Orientation::Positive => {
                    n_over += 1;
                }
                Orientation::Negative => {
                    n_under += 1;
                }
                _ => {}
            }
        }
        if n_over == 0 || n_under == 0 {
            // panic!("don't find plane to split cell");
            return;
        }

        for &eid in &cell_edges {
            let eid = eid.into();
            let e_verts = self.mesh.e_vertices(eid);
            if sign_reversed(
                *self.vert_orientations[tid].get(&e_verts[0]).unwrap(),
                *self.vert_orientations[tid].get(&e_verts[1]).unwrap(),
            ) {
                let new_vid = self.split_edge(eid, &tri, &bump);
                self.vert_orientations[tid].insert(new_vid, Orientation::Zero);
            }
        }

        let n_cell_faces = self.cell_data[cid].faces.len();
        let mut pos_faces = Vec::new();
        let mut neg_faces = Vec::new();
        let mut zero_ori_halfedges = Vec::new_in(bump);
        for i in 0..n_cell_faces {
            let fid = self.cell_data[cid].faces[i];
            let orientaion = &mut self.vert_orientations[tid];
            let coplanar = split_face_verts(&self.mesh, fid, orientaion, bump);
            match coplanar {
                Ok([va, vb]) => {
                    let new_hid = self.mesh.split_face(fid, va, vb, bump);
                    let new_fid = self.mesh.he_face(new_hid);

                    let mut parent = tri.to_vec();

                    parent.extend(&self.face_data[fid].plane);
                    self.edge_data.push(BSPEdgeData::new(parent));
                    let new_is_pos = orientaion
                        .get(&self.mesh.he_tip_vertex(self.mesh.he_next(new_hid)))
                        .unwrap()
                        == &Orientation::Positive;
                    if new_is_pos {
                        pos_faces.push(new_fid);
                        neg_faces.push(fid);
                    } else {
                        pos_faces.push(fid);
                        neg_faces.push(new_fid);
                    }

                    if new_is_pos ^ (self.face_data[fid].cells[0] == cid) {
                        zero_ori_halfedges.push(new_hid);
                    } else {
                        zero_ori_halfedges.push(self.mesh.he_sibling(new_hid));
                    }

                    let old_is_pos = !new_is_pos;

                    for &nei_cid in &self.face_data[fid].cells {
                        if nei_cid != INVALID_IND && nei_cid != cid {
                            self.cell_data[nei_cid].faces.push(new_fid);
                        }
                    }

                    let mut new_face_triangles = Vec::new();
                    self.face_data[fid].triangles.retain(|&face_tid| {
                        let face_tri = {
                            let start = face_tid * 3;
                            &self.constraints[start..(start + 3)]
                        };
                        verts_orient_wrt_plane(
                            &self.points[tri[0]],
                            &self.points[tri[1]],
                            &self.points[tri[2]],
                            face_tri,
                            &self.points,
                            orientaion,
                            bump,
                        );
                        let (mut has_pos, mut has_neg) = (false, false);
                        for vid in face_tri {
                            match orientaion.get(vid).unwrap() {
                                Orientation::Positive => {
                                    has_pos = true;
                                }
                                Orientation::Negative => {
                                    has_neg = true;
                                }
                                _ => {}
                            }
                        }
                        if new_is_pos {
                            if has_pos {
                                new_face_triangles.push(face_tid);
                            }
                        } else {
                            if has_neg {
                                new_face_triangles.push(face_tid);
                            }
                        }
                        if old_is_pos {
                            has_pos
                        } else {
                            has_neg
                        }
                    });

                    let old_face_data = &self.face_data[fid];
                    self.face_data.push(BSPFaceData {
                        cells: old_face_data.cells,
                        triangles: new_face_triangles,
                        plane: old_face_data.plane,
                    });
                    self.edge_visits.push(false);
                }
                Err(coplanar) => match coplanar {
                    Ok(halfedges) => {
                        let vid = self
                            .mesh
                            .he_tip_vertex(self.mesh.he_next(*halfedges.last().unwrap()));
                        // it is impossible orientation is zero
                        let is_pos = *orientaion.get(&vid).unwrap() == Orientation::Positive;
                        if is_pos {
                            pos_faces.push(fid);
                        } else {
                            neg_faces.push(fid);
                        }
                        if is_pos {
                            if self.face_data[fid].cells[0] == cid {
                                // every edge always has two opposite halfedge
                                zero_ori_halfedges.extend(
                                    halfedges.into_iter().map(|hid| self.mesh.he_twin(hid)),
                                );
                            } else {
                                zero_ori_halfedges.extend(halfedges);
                            }
                        }
                    }
                    Err(has_pos) => {
                        if has_pos {
                            pos_faces.push(fid);
                        } else {
                            neg_faces.push(fid);
                        }
                    }
                },
            }
        }

        // add separating face
        {
            let face_halfedges = make_loop(&self.mesh, zero_ori_halfedges);
            let new_fid = self.mesh.add_face_by_halfedges(&face_halfedges, bump);
            let new_cid = self.cell_data.len();
            self.face_data.push(BSPFaceData::new(
                new_cid,
                cid,
                coplanar_triangles,
                [tri[0], tri[1], tri[2]],
            ));

            neg_faces.push(new_fid);
            let [pos_inner_triangles, neg_inner_triangles] =
                self.separate_cell_triangles(cid, tid, bump);
            self.cell_data[cid].faces = neg_faces;
            self.cell_data[cid].inner_triangles = neg_inner_triangles;

            for &fid in &pos_faces {
                for f_cid in &mut self.face_data[fid].cells {
                    if *f_cid == cid {
                        *f_cid = new_cid;
                    }
                }
            }

            pos_faces.push(new_fid);
            self.cell_data
                .push(BSPCellData::new(pos_faces, pos_inner_triangles));
        }

        // if let Err(err) = validate_mesh_connectivity(&self.mesh) {
        //     panic!("the err is {}", err);
        // }

        self.split_duration += start.elapsed();
    }

    pub fn complex_parition(&mut self, tri_in_shell: &[usize]) {
        for cid in 0..self.cell_data.len() {
            self.validate_cell(cid);
        }

        let mut explicit_points = vec![0.0; self.points.len() * 3];
        for (p, data) in self.points.iter().zip(explicit_points.chunks_mut(3)) {
            match p {
                Point3D::Explicit(p) => {
                    data[0] = p.data[0];
                    data[1] = p.data[1];
                    data[2] = p.data[2];
                }
                Point3D::LPI(p) => {
                    p.to_explicit(data);
                }
                Point3D::TPI(p) => {
                    p.to_explicit(data);
                }
            }
        }
        let mut bump = Bump::new();
        let face_areas = {
            let mut areas = Vec::from_iter(self.mesh.faces().map(|fid| {
                bump.reset();
                let mut face_verts = Vec::new_in(&bump);
                face_verts.extend(
                    self.mesh
                        .face(fid)
                        .halfedges()
                        .map(|hid| self.mesh.he_vertex(hid)),
                );
                face_area(&explicit_points, &face_verts, &bump)
            }));
            let area_sum = areas.iter().sum::<f64>();
            for area in &mut areas {
                *area /= area_sum;
            }
            areas
        };
        let n_shells = *tri_in_shell.iter().max().unwrap() + 1;
        let mut cell_costs_external = vec![vec![0.0; self.cell_data.len() + 1]; n_shells];
        let mut cell_costs_internal = vec![vec![0.0; self.cell_data.len() + 1]; n_shells];
        let mut is_black = vec![vec![false; self.mesh.n_faces()]; n_shells];
        for fid in self.mesh.faces() {
            let face_data = &self.face_data[fid];
            let face_triangles = &face_data.triangles;
            if face_triangles.is_empty() {
                continue;
            }

            bump.reset();
            let first_tid = face_triangles[0];
            // uncoplanar vertex
            let vert = [find_uncoplanar_verts(
                triangle(first_tid, &self.constraints),
                face_data.cells[0],
                &mut self.vert_orientations[first_tid],
                &self.points,
                &self.cell_data[face_data.cells[0]].faces,
                &self.mesh,
                &bump,
            )];

            let face_ori = {
                let tri = &face_data.plane;
                // vert must be not coplanar with tri
                orient3d(
                    &self.points[tri[0]],
                    &self.points[tri[1]],
                    &self.points[tri[2]],
                    &self.points[vert[0]],
                    &bump,
                )
            };
            let mut tri_ori_map = HashMap::<usize, i8, _, _>::new_in(&bump);
            for &tid in face_triangles {
                let tri = triangle(tid, &self.constraints);
                verts_orient_wrt_plane(
                    &self.points[tri[0]],
                    &self.points[tri[1]],
                    &self.points[tri[2]],
                    &vert,
                    &self.points,
                    &mut self.vert_orientations[tid],
                    &bump,
                );
                let ori = *self.vert_orientations[tid].get(&vert[0]).unwrap();
                let shell = tri_in_shell[tid];
                let val = tri_ori_map.entry(shell).or_insert(0);
                if ori == face_ori {
                    *val += 1;
                } else {
                    *val -= 1;
                }
            }

            for (shell_id, cnt) in tri_ori_map {
                if cnt > 0 {
                    cell_costs_internal[shell_id][face_data.cells[0]] += face_areas[fid];
                    let second_cid = face_data.cells[1];
                    if second_cid != INVALID_IND {
                        cell_costs_external[shell_id][second_cid] += face_areas[fid];
                    }
                    is_black[shell_id][fid] = true;
                } else if cnt < 0 {
                    cell_costs_external[shell_id][face_data.cells[0]] += face_areas[fid];
                    let second_cid = face_data.cells[1];
                    if second_cid != INVALID_IND {
                        cell_costs_internal[shell_id][second_cid] += face_areas[fid];
                    }
                    is_black[shell_id][fid] = true;
                }
            }
        }

        let mut g = G::new();
        g.external = cell_costs_external[0].clone();
        g.internal = cell_costs_internal[0].clone();
        *g.internal.last_mut().unwrap() = 1.0;

        let mut graphs = Vec::from_iter(
            cell_costs_external
                .into_iter()
                .zip(cell_costs_internal)
                .map(|(external, mut internal)| {
                    internal[self.cell_data.len()] = 1.0;
                    GraphCut::new(&external, &internal)
                }),
        );
        for fid in self.mesh.faces() {
            if self.face_data[fid].triangles.is_empty() {
                let [c1, mut c2] = self.face_data[fid].cells;
                if c2 == INVALID_IND {
                    c2 = self.cell_data.len();
                }

                for shell_id in 0..n_shells {
                    if !is_black[shell_id][fid] {
                        graphs[shell_id].add_edge(c1, c2, face_areas[fid], face_areas[fid]);
                        g.edges.push(Edge::new(c1, c2, face_areas[fid]));
                    }
                }
            }
        }
        let gj = serde_json::to_string(&g).unwrap();
        std::fs::write("123.json", gj).unwrap();
        let mut cell_kept = vec![false; self.cell_data.len() + 1];
        for graph in &mut graphs {
            graph.max_flow();
            for cid in 0..self.cell_data.len() {
                cell_kept[cid] |= !graph.is_sink[cid];
            }
        }

        let mut kept_faces = Vec::new();
        for fid in self.mesh.faces() {
            let [c1, mut c2] = self.face_data[fid].cells;
            if c2 == INVALID_IND {
                c2 = self.cell_data.len();
            }
            if cell_kept[c1] ^ cell_kept[c2] {
                let verts = if cell_kept[c1] {
                    let mut verts = Vec::from_iter(
                        self.mesh
                            .face(fid)
                            .halfedges()
                            .map(|hid| self.mesh.he_vertex(hid).0),
                    );
                    verts.reverse();
                    verts
                } else {
                    Vec::from_iter(
                        self.mesh
                            .face(fid)
                            .halfedges()
                            .map(|hid| self.mesh.he_vertex(hid).0),
                    )
                };
                kept_faces.push(verts);
            }
        }

        let mut txt = "".to_owned();
        for p in explicit_points.chunks(3) {
            txt.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
        }
        for face in &kept_faces {
            txt.push('f');
            for &vid in face {
                txt.push_str(&format!(" {}", vid + 1));
            }
            txt.push('\n');
        }
        std::fs::write("123.obj", txt).unwrap();
    }

    #[inline(always)]
    fn triangle(&self, tid: usize) -> &[VertexId] {
        let start = tid * 3;
        &self.constraints[start..(start + 3)]
    }

    #[inline(always)]
    fn is_virtual(&self, tid: usize) -> bool {
        tid >= self.n_ori_triangles
    }

    fn separate_out_coplanar_triangles<A: Allocator + Copy>(
        &mut self,
        pivot_tid: usize,
        cid: usize,
        bump: A,
    ) -> Vec<usize> {
        let mut plane_pts = Vec::with_capacity_in(3, bump);
        plane_pts.extend(
            self.triangle(pivot_tid)
                .iter()
                .map(|&vid| &self.points[vid]),
        );
        let mut coplanar_triangles = Vec::new();
        self.cell_data[cid].inner_triangles.retain(|&tid| {
            let start = tid * 3;
            let tri = &self.constraints[start..(start + 3)];
            verts_orient_wrt_plane(
                plane_pts[0],
                plane_pts[1],
                plane_pts[2],
                tri,
                &self.points,
                &mut self.vert_orientations[pivot_tid],
                bump,
            );
            let ori = &self.vert_orientations[pivot_tid];
            let coplanar = *ori.get(&tri[0]).unwrap() == Orientation::Zero
                && *ori.get(&tri[1]).unwrap() == Orientation::Zero
                && *ori.get(&tri[2]).unwrap() == Orientation::Zero;
            if coplanar && tid < self.n_ori_triangles {
                coplanar_triangles.push(tid);
            }
            !coplanar
        });
        coplanar_triangles
    }

    #[inline]
    fn cell_verts_and_edges<A: Allocator + Copy>(
        &mut self,
        cid: usize,
        bump: A,
    ) -> (Vec<VertexId, A>, Vec<EdgeId, A>) {
        let mesh = &self.mesh;
        let mut verts = Vec::new_in(bump);
        let mut edges = Vec::new_in(bump);
        for &fid in &self.cell_data[cid].faces {
            for hid in mesh.face(fid.into()).halfedges() {
                let vid = mesh.he_vertex(hid);
                if !self.vert_visits[vid] {
                    self.vert_visits[vid] = true;
                    verts.push(vid);
                }
                let eid = mesh.he_edge(hid);
                if !self.edge_visits[eid] {
                    self.edge_visits[eid] = true;
                    edges.push(eid);
                }
            }
        }
        for &vid in &verts {
            self.vert_visits[vid] = false;
        }
        for &eid in &edges {
            self.edge_visits[eid] = false;
        }
        (verts, edges)
    }

    #[inline]
    fn split_edge<A: Allocator + Copy>(
        &mut self,
        eid: EdgeId,
        tri: &[VertexId],
        bump: A,
    ) -> VertexId {
        let vid = self.mesh.split_edge(eid, bump);
        self.edge_visits.push(false);
        let p0 = self.points[tri[0]].explicit().unwrap().clone();
        let p1 = self.points[tri[1]].explicit().unwrap().clone();
        let p2 = self.points[tri[2]].explicit().unwrap().clone();
        let ori_edge_parents = self.edge_data[eid].parents.clone();
        if ori_edge_parents.len() > 2 {
            self.points.push(three_planes_intersection(
                [&ori_edge_parents[..3], &ori_edge_parents[3..], tri],
                &self.points,
                bump,
            ));
        } else {
            self.points.push(Point3D::LPI(ImplicitPointLPI::new(
                self.points[ori_edge_parents[0]].explicit().unwrap().clone(),
                self.points[ori_edge_parents[1]].explicit().unwrap().clone(),
                p0,
                p1,
                p2,
            )));
        }
        self.vert_visits.push(false);

        self.edge_visits.push(false);
        self.edge_data.push(BSPEdgeData::new(ori_edge_parents));
        vid
    }

    fn separate_cell_triangles<A: Allocator + Copy>(
        &mut self,
        cid: usize,
        pivot_tid: usize,
        bump: A,
    ) -> [Vec<usize>; 2] {
        let cell_triangles = &self.cell_data[cid].inner_triangles;
        let pivot_tri = self.triangle(pivot_tid);
        let pa = &self.points[pivot_tri[0]];
        let pb = &self.points[pivot_tri[1]];
        let pc = &self.points[pivot_tri[2]];
        let mut pos_triangles = Vec::new();
        let mut neg_triangles = Vec::new();
        for &tid in cell_triangles {
            let tri = triangle(tid, &self.constraints);
            let orientation = &mut self.vert_orientations[pivot_tid];
            verts_orient_wrt_plane(pa, pb, pc, tri, &self.points, orientation, bump);
            let (mut has_pos, mut has_neg) = (false, false);
            for vid in tri {
                match orientation.get(vid).unwrap() {
                    Orientation::Positive => {
                        has_pos = true;
                    }
                    Orientation::Negative => {
                        has_neg = true;
                    }
                    _ => {}
                }
            }
            if has_pos {
                pos_triangles.push(tid);
            }
            if has_neg {
                neg_triangles.push(tid);
            }
        }
        [pos_triangles, neg_triangles]
    }

    pub fn validate_cell(&mut self, cid: usize) {
        let bump = Bump::new();
        let (cell_verts, _) = self.cell_verts_and_edges(cid, &bump);
        for &fid in &self.cell_data[cid].faces {
            let mut set = hashbrown::HashSet::new();
            set.extend(
                self.mesh
                    .face(fid)
                    .halfedges()
                    .map(|hid| self.mesh.he_vertex(hid)),
            );
            let face_tri = Vec::from_iter(self.face_data[fid].plane.map(|vid| &self.points[vid]));
            let reversed = self.face_data[fid].cells[0] != cid;
            for &vid in &cell_verts {
                if set.contains(&vid) {
                    continue;
                }
                let ori = orient3d(
                    face_tri[0],
                    face_tri[1],
                    face_tri[2],
                    &self.points[vid],
                    &bump,
                );
                if ori == Orientation::Zero {
                    continue;
                }
                if !((ori == Orientation::Positive) ^ reversed) {
                    panic!("ori reversed");
                }
            }
        }
    }
}

fn triangle(tid: usize, triangles: &[VertexId]) -> &[VertexId] {
    let start = tid * 3;
    &triangles[start..(start + 3)]
}

fn remove_ghost_tets<'b>(mesh: &TetMesh) -> Vec<usize> {
    let tets = &mesh.tets;
    let mut idx = 0;
    let mut new_orders = vec![INVALID_IND; tets.len()];
    for tid in 0..tets.len() {
        if !mesh.is_hull_tet(tid) {
            new_orders[tid] = idx;
            idx += 1;
        }
    }
    new_orders
}

fn insert_coplanar_triangles(src: &[usize], dest: &mut Vec<usize>, n_constraints: usize) {
    for &idx in src {
        if idx < n_constraints {
            if let Err(pos) = dest.binary_search(&idx) {
                dest.insert(pos, idx);
            }
        }
    }
}

#[inline]
fn verts_orient_wrt_plane<A: Allocator + Copy>(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    verts: &[VertexId],
    points: &[Point3D],
    vert_orientations: &mut HashMap<VertexId, Orientation>,
    bump: A,
) {
    let p0 = pa.explicit().unwrap();
    let p1 = pb.explicit().unwrap();
    let p2 = pc.explicit().unwrap();
    for &vid in verts {
        if let hashbrown::hash_map::Entry::Vacant(ori) = vert_orientations.entry(vid) {
            if is_point_built_from_plane(&points[vid], p0, p1, p2) {
                ori.insert(Orientation::Zero);
            } else {
                ori.insert(orient3d(pa, pb, pc, &points[vid], bump));
            }
        }
    }
}

#[inline(always)]
fn is_point_built_from_plane(
    p: &Point3D,
    pa: &ExplicitPoint3D,
    pb: &ExplicitPoint3D,
    pc: &ExplicitPoint3D,
) -> bool {
    match p {
        Point3D::Explicit(p) => {
            if p == pa || p == pb || p == pc {
                return true;
            }
        }
        Point3D::LPI(lpi) => {
            if (((pa == &lpi.p) && (pb == &lpi.q)) || ((pb == &lpi.p) && (pa == &lpi.q)))
                || (((pb == &lpi.p) && (pc == &lpi.q)) || ((pc == &lpi.p) && (pb == &lpi.q)))
                || (((pc == &lpi.p) && (pa == &lpi.q)) || ((pa == &lpi.p) && (pc == &lpi.q)))
                || ((pa == &lpi.r) && (pb == &lpi.s) && (pc == &lpi.t))
            {
                return true;
            }
        }
        Point3D::TPI(tpi) => {
            if ((pa == &tpi.v1) && (pb == &tpi.v2) && (pc == &tpi.v3))
                || ((pa == &tpi.w1) && (pb == &tpi.w2) && (pc == &tpi.w3))
                || ((pa == &tpi.u1) && (pb == &tpi.u2) && (pc == &tpi.u3))
            {
                return true;
            }
        }
    }
    return false;
}

#[inline]
fn three_planes_intersection<A: Allocator + Copy>(
    planes: [&[VertexId]; 3],
    points: &[Point3D],
    bump: A,
) -> Point3D {
    for (&tri, &tri1, &tri2) in planes.iter().circular_tuple_windows() {
        let mut common = Vec::new_in(bump);
        for va in tri1 {
            if tri2.contains(va) {
                common.push(*va);
            }
            if common.len() > 1 {
                break;
            }
        }
        if common.len() > 1 {
            return Point3D::LPI(ImplicitPointLPI::new(
                points[common[0]].explicit().unwrap().clone(),
                points[common[1]].explicit().unwrap().clone(),
                points[tri[0]].explicit().unwrap().clone(),
                points[tri[1]].explicit().unwrap().clone(),
                points[tri[2]].explicit().unwrap().clone(),
            ));
        }
    }
    let [t0, t1, t2] = planes;
    Point3D::TPI(ImplicitPointTPI::new(
        points[t0[0]].explicit().unwrap().clone(),
        points[t0[1]].explicit().unwrap().clone(),
        points[t0[2]].explicit().unwrap().clone(),
        points[t1[0]].explicit().unwrap().clone(),
        points[t1[1]].explicit().unwrap().clone(),
        points[t1[2]].explicit().unwrap().clone(),
        points[t2[0]].explicit().unwrap().clone(),
        points[t2[1]].explicit().unwrap().clone(),
        points[t2[2]].explicit().unwrap().clone(),
    ))
}

fn split_face_verts<A: Allocator + Copy>(
    mesh: &SurfaceMesh,
    fid: FaceId,
    vert_orientations: &HashMap<VertexId, Orientation>,
    bump: A,
) -> Result<[VertexId; 2], Result<Vec<HalfedgeId, A>, bool>> {
    let (mut first, mut second) = (None, None);
    let (mut has_pos, mut has_neg) = (false, false);
    let mut zero_ori_candidates = Vec::new_in(bump);
    for hid in mesh.face(fid).halfedges() {
        let vid = mesh.he_vertex(hid);
        match vert_orientations.get(&vid).unwrap() {
            Orientation::Positive => {
                has_pos = true;
                if has_neg && let Some(first) = first && let Some(second) = second {
                    return Ok([first, second]);
                }
            }
            Orientation::Negative => {
                has_neg = true;
                if has_pos && let Some(first) = first && let Some(second) = second {
                    return Ok([first, second]);
                }
            }
            Orientation::Zero => {
                zero_ori_candidates.push(hid);
                if first.is_none() {
                    first = Some(vid);
                } else if second.is_none() {
                    if has_pos && has_neg {
                        return Ok([first.unwrap(), vid]);
                    } else {
                        second = Some(vid);
                    }
                }
            }
            Orientation::Undefined => {}
        }
    }

    zero_ori_candidates.retain(|&hid| {
        vert_orientations.get(&mesh.he_tip_vertex(hid)).unwrap() == &Orientation::Zero
    });

    if zero_ori_candidates.is_empty() {
        return Err(Err(has_pos));
    } else {
        if zero_ori_candidates.len() == 1 {
            return Err(Ok(zero_ori_candidates));
        } else {
            let pos = zero_ori_candidates
                .windows(2)
                .position(|pair| mesh.he_tip_vertex(pair[0]) != mesh.he_vertex(pair[1]));
            if let Some(pos) = pos {
                zero_ori_candidates.rotate_left(pos);
            }
            return Err(Ok(zero_ori_candidates));
        }
    }
}

fn make_loop<A: Allocator + Copy>(
    mesh: &SurfaceMesh,
    halfedges: Vec<HalfedgeId, A>,
) -> Vec<HalfedgeId, A> {
    let bump = *halfedges.allocator();
    let mut map = HashMap::new_in(bump);
    map.extend(halfedges.iter().map(|&hid| (mesh.he_vertex(hid), hid)));

    let first = halfedges[0];
    let mut result = Vec::new_in(bump);
    result.push(first);
    let mut curr = first;
    loop {
        let next = *map.get(&mesh.he_tip_vertex(curr)).unwrap();
        if next == first {
            break;
        }
        result.push(next);
        curr = next;
    }
    result
}

#[inline(always)]
fn find_uncoplanar_verts<A: Allocator + Copy>(
    tri: &[VertexId],
    cid: usize,
    orientations: &mut HashMap<VertexId, Orientation>,
    points: &[Point3D],
    faces: &[FaceId],
    mesh: &SurfaceMesh,
    bump: A,
) -> VertexId {
    let pa = &points[tri[0]];
    let pb = &points[tri[1]];
    let pc = &points[tri[2]];
    for &fid in faces {
        for vid in mesh.face(fid).halfedges().map(|hid| mesh.he_vertex(hid)) {
            verts_orient_wrt_plane(pa, pb, pc, &[vid], &points, orientations, bump);
            if *orientations.get(&vid).unwrap() != Orientation::Zero {
                return vid;
            }
        }
    }
    panic!("no uncoplanar vertex in cell");
}

#[inline(always)]
fn face_area<A: Allocator + Copy>(points: &[f64], verts: &[VertexId], bump: A) -> f64 {
    let pa = point(points, verts[0].0);
    let mut v1 = Vec::with_capacity_in(3, bump);
    v1.resize(3, 0.0);
    let mut v2 = Vec::with_capacity_in(3, bump);
    v2.resize(3, 0.0);
    let mut n = Vec::with_capacity_in(3, bump);
    n.resize(3, 0.0);
    let mut area = 0.0;
    for two_verts in verts[1..].windows(2) {
        let pb = point(points, two_verts[0].0);
        sub(pb, pa, &mut v1);
        let pc = point(points, two_verts[1].0);
        sub(pc, pa, &mut v2);
        cross(&v1, &v2, &mut n);
        area += norm(&n);
    }
    return area;
}
