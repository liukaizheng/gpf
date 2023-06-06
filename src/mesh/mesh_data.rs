use std::{
    cell::RefCell,
    rc::{Rc, Weak},
};

use super::{ElementId, HalfedgeId, Mesh};
use bumpalo::{collections::Vec, Bump};

pub trait MeshData {
    type Id: ElementId;
    fn expand(&mut self, len: usize);
    fn len(&self) -> usize;
}

pub struct HalfedgeData<'b, T: 'b + Clone, M: Mesh<'b>> {
    default_val: T,
    data: Vec<'b, T>,
    mesh: Weak<RefCell<M>>,
}

impl<'b, T: 'b + Clone, M: Mesh<'b>> HalfedgeData<'b, T, M> {
    #[inline]
    pub fn new(mesh: Weak<RefCell<M>>, bump: &'b Bump, default_val: T) -> Rc<RefCell<Self>> {
        let n_halfedges = mesh.upgrade().map_or(0, |m| m.borrow().n_halfedges());
        let default_val_clone = default_val.clone();
        let data = Rc::new(RefCell::new(Self {
            default_val,
            data: bumpalo::vec![in bump; default_val_clone; n_halfedges],
            mesh: mesh.clone(),
        }));
        if let Some(mesh) = mesh.upgrade() {
            mesh.borrow_mut().add_halfedges_data(Rc::downgrade(&data));
        }
        data
    }
}

impl<'b, T: 'b + Clone, M: Mesh<'b>> HalfedgeData<'b, T, M> {
    fn drop(&mut self) {
        if let Some(mesh) = self.mesh.upgrade() {
            mesh.borrow_mut().remove_halfedges_data(self);
        }
    }
}

impl<'b, T: 'b + Clone, M: Mesh<'b>> MeshData for HalfedgeData<'b, T, M> {
    type Id = HalfedgeId;

    #[inline(always)]
    fn expand(&mut self, len: usize) {
        // self.data.resize(len, Default::default());
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.data.len()
    }
}
