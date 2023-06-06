use std::collections::HashMap;

use bumpalo::collections::Vec;

use crate::{
    predicates::{ExplicitPoint3D, Point3D},
    triangle::TetMesh,
    INVALID_IND,
};

pub(crate) struct BSPComplex<'b> {
    points: Vec<'b, Point3D<'b>>,
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
        Self { points }
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
