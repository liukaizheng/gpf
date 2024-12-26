use std::alloc::Allocator;

use super::ElementId;

pub(super) struct ElementCache<E: ElementId + Copy, A: Allocator + Copy> {
    next_arr: Vec<E, A>,
    curr: E,
}

impl<E: ElementId + Copy, A: Allocator + Copy> ElementCache<E, A> {
    pub(super) fn new(n_elements: usize, alloc: A) -> Self {
        let mut next_arr = Vec::with_capacity_in(n_elements, alloc);
        next_arr.resize(n_elements, E::default());
        Self {
            next_arr,
            curr: E::default(),
        }
    }

    #[inline]
    pub(super) fn push(&mut self, ele: E) {
        self.next_arr[ele.index()] = self.curr;
        self.curr = ele;
    }

    #[inline]
    pub(super) fn pop(&mut self) -> E {
        let ele = self.curr;
        if ele.valid() {
            self.curr = self.next_arr[ele.index()];
        }
        ele
    }

    #[inline]
    pub(super) fn reserve(&mut self, n_elements: usize) {
        self.next_arr.resize(n_elements, E::default());
    }
}
