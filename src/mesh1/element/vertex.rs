use std::ops::Deref;

use crate::INVALID_IND;

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(usize);

impl Default for VertexId {
    #[inline]
    fn default() -> Self {
        VertexId(INVALID_IND)
    }
}

impl From<usize> for VertexId {
    #[inline]
    fn from(index: usize) -> Self {
        VertexId(index)
    }
}

impl ElementId for VertexId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for VertexId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
