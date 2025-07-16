use std::ops::Deref;

use crate::INVALID_IND;

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfedgeId(pub usize);

impl Default for HalfedgeId {
    #[inline]
    fn default() -> Self {
        HalfedgeId(INVALID_IND)
    }
}

impl From<usize> for HalfedgeId {
    #[inline]
    fn from(index: usize) -> Self {
        HalfedgeId(index)
    }
}

impl ElementId for HalfedgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for HalfedgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
