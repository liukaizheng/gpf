use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{iter_next, Element, ElementId};
use crate::{element_id, mesh::ManifoldMesh, INVALID_IND};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoopId(usize);

element_id! {struct BoundaryLoopId}

pub struct BoundaryLoopIter<'a> {
    id: BoundaryLoopId,
    mesh: &'a ManifoldMesh,
}

impl<'a> BoundaryLoopIter<'a> {
    pub fn new(id: BoundaryLoopId, mesh: &'a ManifoldMesh) -> Self {
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

impl<'a> Deref for BoundaryLoopIter<'a> {
    type Target = BoundaryLoopId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a> Element for BoundaryLoopIter<'a> {
    type Id = BoundaryLoopId;
    type M = ManifoldMesh;

    fn id(&self) -> BoundaryLoopId {
        self.id
    }

    fn mesh(&self) -> &Self::M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.boundary_loop_is_valid(self.id)
    }

    fn next(&mut self) {
        self.id.0 += 1;
    }

    fn is_end(&self) -> bool {
        self.id.0 == self.capacity()
    }
}

impl<'a> Iterator for BoundaryLoopIter<'a> {
    type Item = BoundaryLoopId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
