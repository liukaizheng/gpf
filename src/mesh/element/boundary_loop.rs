use std::ops::{Add, Deref, DerefMut, Index, IndexMut, Mul};

use super::{ElementId, ElementIndex};
use crate::{element_id, INVALID_IND};

use std::alloc::Allocator;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryLoopId(pub usize);

element_id!(struct BoundaryLoopId);
