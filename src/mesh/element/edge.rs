use std::ops::{Deref, Index, IndexMut};

use super::{iter_next, Element, HalfedgeIter};
use crate::mesh::Mesh;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

impl From<usize> for EdgeId {
    fn from(id: usize) -> Self {
        Self(id)
    }
}

impl<T> Index<EdgeId> for Vec<T> {
    type Output = T;
    fn index(&self, index: EdgeId) -> &Self::Output {
        &self[index.0]
    }
}

impl<T> IndexMut<EdgeId> for Vec<T> {
    fn index_mut(&mut self, index: EdgeId) -> &mut Self::Output {
        &mut self[index.0]
    }
}

pub struct EdgeIter<'a, M: Mesh> {
    id: EdgeId,
    mesh: &'a M,
}

impl<'a, M: Mesh> EdgeIter<'a, M> {
    pub fn new(id: EdgeId, mesh: &'a M) -> Self {
        Self { id, mesh }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.mesh.n_halfedges()
    }

    fn capacity(&self) -> usize {
        self.mesh.n_halfedges_capacity()
    }
}

impl<'a, M: Mesh> Deref for EdgeIter<'a, M> {
    type Target = EdgeId;

    fn deref(&self) -> &Self::Target {
        &self.id
    }
}

impl<'a, M: Mesh> Element for EdgeIter<'a, M> {
    type Item = EdgeId;
    type M = M;

    fn id(&self) -> EdgeId {
        self.id
    }

    fn mesh(&self) -> &M {
        self.mesh
    }

    fn valid(&self) -> bool {
        self.mesh.edge_is_valid(self.id)
    }

    fn next(&mut self) -> bool {
        self.id.0 += 1;
        self.id.0 < self.capacity()
    }
}

pub trait Edge: Element {
    fn halfedge(&self) -> HalfedgeIter<'_, Self::M>;
}

impl<'a, M: Mesh> Edge for EdgeIter<'a, M> {
    fn halfedge(&self) -> HalfedgeIter<'_, Self::M> {
        HalfedgeIter::new(self.mesh.e_halfedge(self.id), self.mesh)
    }
}

impl<'a, M: Mesh> Iterator for EdgeIter<'a, M> {
    type Item = EdgeId;

    fn next(&mut self) -> Option<Self::Item> {
        iter_next(self)
    }
}
