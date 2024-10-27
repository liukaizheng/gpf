use std::ops::{Deref, DerefMut, Index, IndexMut};

use super::{ElementId, ElementIndex};
use crate::{element_id, INVALID_IND};

use std::alloc::Allocator;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoopId(pub usize);

element_id!(struct BoundaryLoopId);
