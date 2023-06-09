use std::{cell::RefCell, rc::Rc};

use gpf::mesh::{HalfedgeData, MeshData, SurfaceMesh};

#[test]
fn build_mesh_data() {
    let mesh = Rc::new(RefCell::new(SurfaceMesh::from(vec![
        vec![0, 1, 2],
        vec![0, 2, 3],
        vec![0, 3, 1],
        vec![0, 4, 5],
        vec![0, 5, 6],
        vec![0, 6, 4],
    ])));
    {
        let halfedge_data = HalfedgeData::new(Rc::downgrade(&mesh), 0.0);
        assert_eq!(halfedge_data.borrow().len(), 18);
        assert_eq!(mesh.borrow().halfedges_data.len(), 1);
    }
    assert_eq!(mesh.borrow_mut().halfedges_data.len(), 0);
}
