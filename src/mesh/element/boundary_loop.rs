use std::ops::{Deref, Index, IndexMut};

use super::{iter_next, Element};
use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoopId(usize);

impl From<usize> for BoundaryLoopId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl<T> Index<BoundaryLoopId> for Vec<T> {
    type Output = T;
    fn index(&self, index: BoundaryLoopId) -> &Self::Output {
        &self[index.0]
    }
}

impl<T> IndexMut<BoundaryLoopId> for Vec<T> {
    fn index_mut(&mut self, index: BoundaryLoopId) -> &mut Self::Output {
        &mut self[index.0]
    }
}

pub struct BoundaryLoopIter<'a, M: Mesh> {
    id: BoundaryLoopId,
    mesh: &'a M,
}

impl<'a, M: Mesh> BoundaryLoopIter<'a, M> {
    pub fn new(id: BoundaryLoopId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_boundary_loops()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_boundary_loop_capacity()
    }
}

impl<'a, M: Mesh> Deref for BoundaryLoopIter<'a, M> {
    type Target = BoundaryLoopId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, M: Mesh> Element for BoundaryLoopIter<'a, M> {
    type Item = BoundaryLoopId;
    type M = M;

    fn id(&self) -> BoundaryLoopId {
        self.id
    }

    fn mesh(&self) -> &M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.boundary_loop_is_valid(self.id)
    }

    fn next(&mut self) -> bool {
        self.id.0 += 1;
        self.id.0 < self.capacity()
    }
}

impl<'a, M: Mesh> Iterator for BoundaryLoopIter<'a, M> {
    type Item = BoundaryLoopId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
