use std::ops::Deref;

use crate::INVALID_IND;

use super::ElementId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

impl Default for FaceId {
    #[inline]
    fn default() -> Self {
        FaceId(INVALID_IND)
    }
}

impl From<usize> for FaceId {
    #[inline]
    fn from(index: usize) -> Self {
        FaceId(index)
    }
}

impl ElementId for FaceId {
    #[inline]
    fn index(&self) -> usize {
        self.0
    }
}

impl Deref for FaceId {
    type Target = usize;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
