use std::{cell::RefCell, rc::Weak};

use super::{ElementId, HalfedgeId, Mesh};

pub trait MeshData {
    type Id: ElementId;
    fn expand(&mut self, len: usize);
    fn len(&self) -> usize;
}

pub struct HalfedgeData<T: Default + Clone, M: Mesh> {
    data: Vec<T>,
    mesh: Weak<RefCell<M>>,
}

impl<T: Default + Clone, M: Mesh> HalfedgeData<T, M> {
    #[inline]
    pub fn new(mesh: Weak<RefCell<M>>) -> Self {
        let n_halfedges = mesh.upgrade().map_or(0, |m| m.borrow().n_halfedges());
        Self {
            data: vec![T::default(); n_halfedges],
            mesh,
        }
    }
}

impl<T: Default + Clone, M: Mesh> Drop for HalfedgeData<T, M> {
    fn drop(&mut self) {
        if let Some(mesh) = self.mesh.upgrade() {
            mesh.borrow_mut().remove_halfedges_data(self);
        }
    }
}

impl<T: Default + Clone, M: Mesh> MeshData for HalfedgeData<T, M> {
    type Id = HalfedgeId;

    #[inline(always)]
    fn expand(&mut self, len: usize) {
        self.data.resize(len, Default::default());
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.data.len()
    }
}
