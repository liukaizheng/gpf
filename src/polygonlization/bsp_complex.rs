use std::{cell::RefCell, rc::Rc};

use bumpalo::collections::Vec;

use crate::{
    mesh::{FaceData, SurfaceMesh, VertexData},
    predicates::{ExplicitPoint3D, Point3D},
    triangle::TetMesh,
    INVALID_IND,
};

use super::conforming_mesh::Constraints;

#[derive(Clone)]
struct BSPFaceData<'b> {
    cells: [usize; 2],
    // coplanar triangles
    triangles: Vec<'b, usize>,
}

impl<'b> BSPFaceData<'b> {
    fn new(c1: usize, c2: usize, triangles: Vec<'b, usize>) -> Self {
        Self {
            cells: [c1, c2],
            triangles,
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
    points: Rc<RefCell<VertexData<'b, Point3D<'b>, SurfaceMesh<'b>>>>,
    mesh: Rc<RefCell<SurfaceMesh<'b>>>,
    face_data: Rc<RefCell<FaceData<'b, BSPFaceData<'b>, SurfaceMesh<'b>>>>,
    constraints: &'a Constraints<'b>,
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
                    insert_coplanar_triangles(
                        &tet_mark[j][nei.tet],
                        &mut face_data[fid].triangles,
                        constraints.n_ori_triangles,
                    );
                    cell_faces.push(fid);
                }
            }
            cell_data.push(BSPCellData::new(cell_faces, tet_mark[4][tid].clone()))
        }

        let mesh = Rc::new(RefCell::new(SurfaceMesh::from(faces)));
        let face_data = FaceData::from_data(
            Rc::downgrade(&mesh),
            face_data,
            BSPFaceData::new(INVALID_IND, INVALID_IND, bumpalo::vec![in bump]),
        );
        let points = VertexData::from_data(
            Rc::downgrade(&mesh),
            points,
            Point3D::Explicit(ExplicitPoint3D {
                data: [0.0, 0.0, 0.0],
            }),
        );

        Self {
            points,
            mesh,
            face_data,
            constraints,
        }
    }

    fn split_cell(cid: usize) {}
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
