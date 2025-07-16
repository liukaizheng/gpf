
use std::ops::Deref;

use crate::INVALID_IND;

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

impl Default for EdgeId {
    #[inline]
    fn default() -> Self {
        EdgeId(INVALID_IND)
    }
}

impl From<usize> for EdgeId {
    #[inline]
    fn from(index: usize) -> Self {
        EdgeId(index)
    }
}

impl ElementId for EdgeId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for EdgeId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
