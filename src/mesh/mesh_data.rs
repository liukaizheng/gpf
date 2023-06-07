use std::{
    cell::RefCell,
    rc::{Rc, Weak},
};

use super::{EdgeId, ElementId, FaceId, HalfedgeId, Mesh, VertexId};
use bumpalo::collections::Vec;

pub trait MeshData {
    type Id: ElementId;
    fn expand(&mut self, len: usize);
    fn len(&self) -> usize;
}

macro_rules! mesh_data {
    (struct $name: ident, $n_elements: ident, $add_ele_data: ident, $remove_ele_data: ident, $id: ty) => {
        pub struct $name<'b, T: 'b + Clone, M: Mesh<'b>> {
            default_val: T,
            data: Vec<'b, T>,
            mesh: Weak<RefCell<M>>,
        }

        impl<'b, T: 'b + Clone, M: Mesh<'b>> $name<'b, T, M> {
            #[inline]
            pub fn new(mesh: Weak<RefCell<M>>, default_val: T) -> Rc<RefCell<Self>> {
                let m = mesh.upgrade().unwrap();
                let n_elements = m.borrow().$n_elements();
                let bump = m.borrow().bump();
                let default_val_clone = default_val.clone();
                let data = Rc::new(RefCell::new(Self {
                    default_val,
                    data: bumpalo::vec![in bump; default_val_clone; n_elements],
                    mesh: mesh.clone(),
                }));
                if let Some(mesh) = mesh.upgrade() {
                    mesh.borrow_mut().$add_ele_data(Rc::downgrade(&data));
                }
                data
            }
        }
        impl<'b, T: 'b + Clone, M: Mesh<'b>> Drop for $name<'b, T, M> {
            fn drop(&mut self) {
                if let Some(mesh) = self.mesh.upgrade() {
                    mesh.borrow_mut().$remove_ele_data(self);
                }
            }
        }

        impl<'b, T: 'b + Clone, M: Mesh<'b>> MeshData for $name<'b, T, M> {
            type Id = $id;

            #[inline(always)]
            fn expand(&mut self, len: usize) {
                self.data.resize(len, self.default_val.clone());
            }

            #[inline(always)]
            fn len(&self) -> usize {
                self.data.len()
            }
        }
    }
}

mesh_data! {struct VertexData, n_vertices_capacity, add_vertices_data, remove_vertices_data, VertexId}
mesh_data! {struct HalfedgeData, n_halfedges_capacity, add_halfedges_data, remove_halfedges_data, HalfedgeId}
mesh_data! {struct EdgeData, n_edges_capacity, add_edges_data, remove_edges_data, EdgeId}
mesh_data! {struct FaceData, n_faces_capacity, add_faces_data, remove_faces_data, FaceId}
