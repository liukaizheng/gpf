use std::{cell::RefCell, collections::HashMap, rc::Rc};

use bumpalo::collections::Vec;

use crate::{
    mesh::{FaceData, SurfaceMesh},
    predicates::{ExplicitPoint3D, Point3D},
    triangle::TetMesh,
    INVALID_IND,
};

pub(crate) struct BSPComplex<'b> {
    points: Vec<'b, Point3D<'b>>,
    mesh: Rc<RefCell<SurfaceMesh<'b>>>,
    face_cells: Rc<RefCell<FaceData<'b, [usize; 2], SurfaceMesh<'b>>>>,
}

#[inline(always)]
fn tet_face_is_new(tid: usize, adj_tid: usize, adj_cid: usize) -> bool {
    adj_tid > tid || adj_cid == INVALID_IND
}

impl<'b> BSPComplex<'b> {
    pub(crate) fn new(tet_mesh: TetMesh<'_, 'b>, tet_mark: [Vec<'b, Vec<'b, usize>>; 5]) -> Self {
        let bump = tet_mesh.tets.bump();
        let points = Vec::from_iter_in(
            tet_mesh
                .points
                .chunks(3)
                .map(|data| Point3D::Explicit(ExplicitPoint3D::from(data))),
            bump,
        );
        let (n_cells, new_tet_orders) = remove_ghost_tets(&tet_mesh);
        let mut face_cells: Vec<'b, [usize; 2]> = Vec::new_in(bump);
        let mut faces: Vec<'b, Vec<'b, usize>> = Vec::new_in(bump);
        for tid in 0..new_tet_orders.len() {
            let cell_idx = new_tet_orders[tid];
            if cell_idx == INVALID_IND {
                continue;
            }
            for nei in &tet_mesh.tets[tid].nei {
                let adj_cell_idx = new_tet_orders[nei.tet];
                if tet_face_is_new(tid, nei.tet, adj_cell_idx) {
                    faces.push(bumpalo::vec![in bump;
                        tet_mesh.org(nei),
                        tet_mesh.dest(nei),
                        tet_mesh.apex(nei)
                    ]);
                    face_cells.push([cell_idx, adj_cell_idx]);
                }
            }
        }
        let mesh = Rc::new(RefCell::new(SurfaceMesh::from(faces)));
        let face_cells =
            FaceData::from_data(Rc::downgrade(&mesh), face_cells, [INVALID_IND, INVALID_IND]);
        Self {
            points,
            mesh,
            face_cells,
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
