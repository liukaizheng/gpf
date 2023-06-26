use std::time::{Duration, Instant};

use bumpalo::Bump;
use hashbrown::HashMap;
use itertools::Itertools;

use crate::{
    mesh::{EdgeId, Face, FaceId, HalfedgeId, Mesh, SurfaceMesh, VertexId},
    predicates::{
        orient3d::orient3d, sign_reversed, ExplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI,
        Orientation, Point3D,
    },
    triangle::TetMesh,
    INVALID_IND,
};

use super::conforming_mesh::Constraints;

#[derive(Clone)]
struct BSPEdgeData {
    parents: Vec<VertexId>,
}

impl BSPEdgeData {
    fn new(parents: Vec<VertexId>) -> Self {
        Self { parents }
    }
}

// #[derive(Clone)]
// enum FaceColor {
// White,
// Black,
// Gray,
// }

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
                    let tri = vec![tet_mesh.org(nei), tet_mesh.dest(nei), tet_mesh.apex(nei)];
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

    pub(crate) fn split_cell(&mut self, cid: usize, bump: &Bump) {
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
        for vid in &cell_verts {
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
            panic!("don't find plane to split cell");
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
        let mut pos_faces = Vec::new_in(bump);
        let mut neg_faces = Vec::new_in(bump);
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
                    let new_is_pos = orientaion.get(&self.mesh.he_tip_vertex(new_hid)).unwrap()
                        == &Orientation::Positive;
                    if new_is_pos {
                        pos_faces.push(new_fid);
                        neg_faces.push(fid);
                    } else {
                        pos_faces.push(fid);
                        neg_faces.push(new_fid);
                    }
                    let old_is_pos = !new_is_pos;

                    for nei_cid in self.face_data[fid].cells {
                        if nei_cid != INVALID_IND {
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
                        let vid = self.mesh.he_tip_vertex(*halfedges.last().unwrap());
                        // it is impossible orientation is zero
                        if orientaion.get(&vid).unwrap() == &Orientation::Positive {
                            pos_faces.push(fid);
                        } else {
                            neg_faces.push(fid);
                        }
                        zero_ori_halfedges.extend(halfedges);
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

        // if let Err(err) = validate_mesh_connectivity(&self.mesh) {
        //     panic!("the err is {}", err);
        // }

        self.split_duration += start.elapsed();
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

    fn separate_out_coplanar_triangles<'b>(
        &mut self,
        pivot_tid: usize,
        cid: usize,
        bump: &'b Bump,
    ) -> Vec<usize, &'b Bump> {
        let mut plane_pts = Vec::with_capacity_in(3, bump);
        plane_pts.extend(
            self.triangle(pivot_tid)
                .iter()
                .map(|&vid| &self.points[vid]),
        );
        let mut coplanar_triangles = Vec::new_in(bump);
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
    fn cell_verts_and_edges<'b>(
        &mut self,
        cid: usize,
        bump: &'b Bump,
    ) -> (Vec<VertexId, &'b Bump>, Vec<EdgeId, &'b Bump>) {
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
    fn split_edge(&mut self, eid: EdgeId, tri: &[VertexId], bump: &Bump) -> VertexId {
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
fn verts_orient_wrt_plane(
    pa: &Point3D,
    pb: &Point3D,
    pc: &Point3D,
    verts: &[VertexId],
    points: &[Point3D],
    vert_orientations: &mut HashMap<VertexId, Orientation>,
    bump: &Bump,
) {
    let p0 = pa.explicit().unwrap();
    let p1 = pb.explicit().unwrap();
    let p2 = pc.explicit().unwrap();
    for &vid in verts {
        if let hashbrown::hash_map::Entry::Vacant(ori) = vert_orientations.entry(vid) {
            if is_point_built_from_plane(&points[vid], p0, p1, p2) {
                ori.insert(Orientation::Zero);
            } else {
                ori.insert(orient3d(&points[vid], pa, pb, pc, bump));
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
fn three_planes_intersection<'b>(
    planes: [&[VertexId]; 3],
    points: &[Point3D],
    bump: &Bump,
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

fn split_face_verts<'b>(
    mesh: &SurfaceMesh,
    fid: FaceId,
    vert_orientations: &HashMap<VertexId, Orientation>,
    bump: &'b Bump,
) -> Result<[VertexId; 2], Result<Vec<HalfedgeId, &'b Bump>, bool>> {
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
