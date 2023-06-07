use std::{cell::RefCell, rc::Rc};

use gpf::mesh::{HalfedgeData, MeshData, SurfaceMesh};

#[test]
fn build_mesh_data() {
    let bump = bumpalo::Bump::new();
    let mesh = Rc::new(RefCell::new(SurfaceMesh::from(bumpalo::vec![in &bump;
        bumpalo::vec![in &bump; 0, 1, 2],
        bumpalo::vec![in &bump; 0, 2, 3],
        bumpalo::vec![in &bump; 0, 3, 1],
        bumpalo::vec![in &bump; 0, 4, 5],
        bumpalo::vec![in &bump; 0, 5, 6],
        bumpalo::vec![in &bump; 0, 6, 4],
    ])));
    {
        let halfedge_data = HalfedgeData::new(Rc::downgrade(&mesh), 0.0);
        let _halfedge_data = HalfedgeData::new(Rc::downgrade(&mesh), 0);
        assert_eq!(halfedge_data.borrow().len(), 18);
        assert_eq!(mesh.borrow().halfedges_data.len(), 2);
    }
    assert_eq!(mesh.borrow_mut().halfedges_data.len(), 0);
}
