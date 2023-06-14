use std::{cell::RefCell, rc::Rc};

use bumpalo::{collections::Vec, Bump};
use itertools::Itertools;

use crate::{
    mesh::{EdgeData, EdgeId, Face, FaceData, Mesh, SurfaceMesh, validate_mesh_connectivity},
    predicates::{
        orient3d::orient3d, sign_reversed, ExplicitPoint3D, ImplicitPointLPI, ImplicitPointTPI,
        Orientation, Point3D,
    },
    triangle::TetMesh,
    INVALID_IND,
};

use super::conforming_mesh::Constraints;

#[derive(Clone)]
struct BSPEdgeData<'b> {
    parents: Vec<'b, usize>,
}

impl<'b> BSPEdgeData<'b> {
    fn new(parents: Vec<'b, usize>) -> Self {
        Self { parents }
    }
}

#[derive(Clone)]
enum FaceColor {
    // White,
    // Black,
    Gray,
}

#[derive(Clone)]
struct BSPFaceData<'b> {
    cells: [usize; 2],
    // coplanar triangles
    triangles: Vec<'b, usize>,
    color: FaceColor,
}

impl<'b> BSPFaceData<'b> {
    fn new(c1: usize, c2: usize, triangles: Vec<'b, usize>) -> Self {
        Self {
            cells: [c1, c2],
            triangles,
            color: FaceColor::Gray,
        }
    }
}

#[derive(Clone)]
struct BSPCellData<'b> {
    faces: Vec<'b, usize>,
    inner_triangles: Vec<'b, usize>,
}

impl<'b> BSPCellData<'b> {
    fn new(faces: Vec<'b, usize>, inner_triangles: Vec<'b, usize>) -> Self {
        Self {
            faces,
            inner_triangles,
        }
    }
}

pub(crate) struct BSPComplex<'a, 'b> {
    points: Vec<'b, Point3D<'b>>,
    mesh: Rc<RefCell<SurfaceMesh<'b>>>,
    edge_data: Rc<RefCell<EdgeData<'b, BSPEdgeData<'b>, SurfaceMesh<'b>>>>,
    face_data: Rc<RefCell<FaceData<'b, BSPFaceData<'b>, SurfaceMesh<'b>>>>,
    cell_data: Vec<'b, BSPCellData<'b>>,
    constraints: &'a Constraints<'b>,

    vert_orientations: Vec<'b, Orientation>,
    vert_visits: Vec<'b, bool>,

    edge_visits: Vec<'b, bool>,
}

#[inline(always)]
fn tet_face_is_new(tid: usize, adj_tid: usize, adj_cid: usize) -> bool {
    adj_tid > tid || adj_cid == INVALID_IND
}

impl<'a, 'b> BSPComplex<'a, 'b> {
    pub(crate) fn new(
        tet_mesh: TetMesh<'_, 'b>,
        constraints: &'a Constraints<'b>,
        tet_mark: [Vec<'b, Vec<'b, usize>>; 5],
    ) -> Self {
        let bump = tet_mesh.tets.bump();
        let points = Vec::from_iter_in(
            tet_mesh
                .points
                .chunks(3)
                .map(|data| Point3D::Explicit(ExplicitPoint3D::from(data))),
            bump,
        );

        let (_, new_tet_orders) = remove_ghost_tets(&tet_mesh);
        let mut face_positions = bumpalo::vec![in bump; INVALID_IND; tet_mesh.tets.len() << 2];
        let mut faces: Vec<'b, Vec<'b, usize>> = Vec::new_in(bump);
        let mut face_data: Vec<'b, BSPFaceData<'b>> = Vec::new_in(bump);
        let mut cell_data: Vec<'b, BSPCellData<'b>> = Vec::new_in(bump);
        for (tid, tet) in tet_mesh.tets.iter().enumerate() {
            let cell_idx = new_tet_orders[tid];
            if cell_idx == INVALID_IND {
                continue;
            }

            let mut cell_faces = Vec::new_in(bump);
            for (i, nei) in tet.nei.iter().enumerate() {
                let adj_cell_idx = new_tet_orders[nei.tet];
                if tet_face_is_new(tid, nei.tet, adj_cell_idx) {
                    faces.push(bumpalo::vec![in bump;
                        tet_mesh.org(nei),
                        tet_mesh.dest(nei),
                        tet_mesh.apex(nei)
                    ]);
                    let mut face = BSPFaceData::new(cell_idx, adj_cell_idx, bumpalo::vec![in bump]);
                    insert_coplanar_triangles(
                        &tet_mark[i][tid],
                        &mut face.triangles,
                        constraints.n_ori_triangles,
                    );
                    let fid = face_data.len();
                    face_positions[(tid << 2) + i] = fid;
                    cell_faces.push(fid);
                    face_data.push(face);
                } else {
                    let j = nei.ver & 3;
                    let fid = face_positions[(nei.tet << 2) + j];
                    let face = &mut face_data[fid];
                    insert_coplanar_triangles(
                        &tet_mark[j][nei.tet],
                        &mut face.triangles,
                        constraints.n_ori_triangles,
                    );
                    cell_faces.push(fid);
                }
            }
            cell_data.push(BSPCellData::new(cell_faces, tet_mark[4][tid].clone()))
        }

        let mesh = SurfaceMesh::from(faces);

        let vert_orientations = bumpalo::vec![in bump; Orientation::Undefined; points.len()];
        let vert_visits = bumpalo::vec![in bump; false; points.len()];

        let edge_visits = bumpalo::vec![in bump; false; mesh.n_edges()];
        let edge_data = Vec::from_iter_in(
            mesh.edges().map(|eid| {
                BSPEdgeData::new(Vec::from_iter_in(
                    mesh.e_vertices(eid).map(|vid| vid.0),
                    bump,
                ))
            }),
            bump,
        );

        let mesh = Rc::new(RefCell::new(mesh));

        let edge_data = EdgeData::from_data(
            Rc::downgrade(&mesh),
            edge_data,
            BSPEdgeData::new(Vec::new_in(bump)),
        );
        let face_data = FaceData::from_data(
            Rc::downgrade(&mesh),
            face_data,
            BSPFaceData::new(INVALID_IND, INVALID_IND, bumpalo::vec![in bump]),
        );
        Self {
            points,
            mesh,
            edge_data,
            face_data,
            cell_data,
            constraints,
            edge_visits,
            vert_orientations,
            vert_visits,
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

    pub(crate) fn split_cell(&mut self, cid: usize) {
        let bump = self.points.bump();
        let tid = *self.cell_data[cid].inner_triangles.last().unwrap();
        let tri = self.constraints.triangle(tid);
        self.cell_data[cid].inner_triangles.pop();
        let mut coplanar_triangles = separate_out_coplanar_triangles(
            tid,
            &mut self.cell_data[cid].inner_triangles,
            &self.points,
            &mut self.vert_orientations,
            &self.constraints,
        );
        if !self.constraints.is_virtual(tid) {
            coplanar_triangles.push(tid);
        }
        let (cell_verts, cell_edges) = self.cell_verts_and_edges(cid);
        verts_orient_wrt_plane(
            &self.points[tri[0]],
            &self.points[tri[1]],
            &self.points[tri[2]],
            &cell_verts,
            &self.points,
            &mut self.vert_orientations,
            bump,
        );

        let (mut n_over, mut n_under) = (0, 0);
        for &vid in &cell_verts {
            match self.vert_orientations[vid] {
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
            unreachable!("don't find plane to split cell");
        }

        for &eid in &cell_edges {
            let eid = eid.into();
            let e_verts = self.mesh.borrow().e_vertices(eid);
            if sign_reversed(
                self.vert_orientations[e_verts[0]],
                self.vert_orientations[e_verts[1]],
            ) {
                self.split_edge(eid, tri);
                let res = validate_mesh_connectivity(&self.mesh.borrow());
                if let Err(str) = res {
                    panic!("failed to split {}", str);
                }
            }
        }

        self.vert_orientations.fill(Orientation::Undefined);
    }

    #[inline]
    fn cell_verts_and_edges(&mut self, cid: usize) -> (Vec<'b, usize>, Vec<'b, EdgeId>) {
        let mesh = self.mesh.borrow();
        let mut verts = Vec::new_in(mesh.bump());
        let mut edges = Vec::new_in(mesh.bump());
        for &fid in &self.cell_data[cid].faces {
            for hid in mesh.face(fid.into()).halfedges() {
                let vid = mesh.he_vertex(hid).0;
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
    fn split_edge(&mut self, eid: EdgeId, tri: &[usize]) {
        self.mesh.borrow_mut().split_edge(eid);
        self.edge_visits.push(false);
        let p0 = self.points[tri[0]].explicit().unwrap().clone();
        let p1 = self.points[tri[1]].explicit().unwrap().clone();
        let p2 = self.points[tri[2]].explicit().unwrap().clone();
        let ori_edge_parents = self.edge_data.borrow().data[eid].parents.clone();
        if ori_edge_parents.len() > 2 {
            self.points.push(three_planes_intersection(
                [&ori_edge_parents[..3], &ori_edge_parents[3..], tri],
                &self.points,
                ori_edge_parents.bump(),
            ));
        } else {
            self.points.push(Point3D::LPI(ImplicitPointLPI::new(
                self.points[ori_edge_parents[0]].explicit().unwrap().clone(),
                self.points[ori_edge_parents[1]].explicit().unwrap().clone(),
                p0,
                p1,
                p2,
                ori_edge_parents.bump(),
            )));
        }
        self.vert_orientations.push(Orientation::Zero);
        self.vert_visits.push(false);

        self.edge_data.borrow_mut().data.last_mut().unwrap().parents = ori_edge_parents;
    }
}

fn remove_ghost_tets<'b>(mesh: &TetMesh<'_, 'b>) -> (usize, Vec<'b, usize>) {
    let tets = &mesh.tets;
    let mut idx = 0;
    let mut new_orders = bumpalo::vec![in tets.bump(); INVALID_IND; tets.len()];
    for tid in 0..tets.len() {
        if !mesh.is_hull_tet(tid) {
            new_orders[tid] = idx;
            idx += 1;
        }
    }
    (idx, new_orders)
}

fn insert_coplanar_triangles<'b>(src: &[usize], dest: &mut Vec<'b, usize>, n_constraints: usize) {
    for &idx in src {
        if idx < n_constraints {
            if let Err(pos) = dest.binary_search(&idx) {
                dest.insert(pos, idx);
            }
        }
    }
}

fn separate_out_coplanar_triangles<'b>(
    pivot_tid: usize,
    triangles: &mut Vec<'b, usize>,
    points: &[Point3D<'b>],
    vert_orientations: &mut [Orientation],
    constraints: &Constraints,
) -> Vec<'b, usize> {
    let bump = triangles.bump();
    let plane_pts = Vec::from_iter_in(
        constraints
            .triangle(pivot_tid)
            .iter()
            .map(|&vid| &points[vid]),
        bump,
    );
    let mut coplanar_triangles: Vec<usize> = Vec::new_in(bump);
    triangles.retain(|&tid| {
        let tri = constraints.triangle(tid);
        verts_orient_wrt_plane(
            plane_pts[0],
            plane_pts[1],
            plane_pts[2],
            tri,
            points,
            vert_orientations,
            bump,
        );
        let coplanar = vert_orientations[tri[0]] == Orientation::Zero
            && vert_orientations[tri[1]] == Orientation::Zero
            && vert_orientations[tri[2]] == Orientation::Zero;
        if coplanar && !constraints.is_virtual(tid) {
            coplanar_triangles.push(tid);
        }
        !coplanar
    });
    coplanar_triangles
}

#[inline(always)]
fn verts_orient_wrt_plane<'b>(
    pa: &Point3D<'b>,
    pb: &Point3D<'b>,
    pc: &Point3D<'b>,
    verts: &[usize],
    points: &[Point3D<'b>],
    vert_orientations: &mut [Orientation],
    bump: &'b Bump,
) {
    let p0 = pa.explicit().unwrap();
    let p1 = pb.explicit().unwrap();
    let p2 = pc.explicit().unwrap();
    for &vid in verts {
        if vert_orientations[vid] != Orientation::Undefined {
            continue;
        }
        if is_point_built_from_plane(&points[vid], p0, p1, p2) {
            vert_orientations[vid] = Orientation::Zero;
            continue;
        }
        vert_orientations[vid] = orient3d(&points[vid], pa, pb, pc, bump);
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
    planes: [&[usize]; 3],
    points: &[Point3D<'b>],
    bump: &'b Bump,
) -> Point3D<'b> {
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
                bump,
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
        bump,
    ))
}
