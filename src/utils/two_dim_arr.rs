use std::{alloc::Allocator, ops::Index};

pub struct TwoDimArr<T, A: Allocator + Copy = std::alloc::Global> {
    pub data: Vec<T, A>,
    pub separators: Vec<usize, A>,
}

impl<T, A: Allocator + Copy> TwoDimArr<T, A> {
    #[inline]
    pub fn new_in(alloc: A) -> Self {
        let mut separators = Vec::new_in(alloc);
        separators.push(0);
        Self {
            data: Vec::new_in(alloc),
            separators,
        }
    }

    #[inline]
    pub fn iter(&self) -> TwoDimArrIter<T, A> {
        TwoDimArrIter { arr: self, idx: 0 }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.separators.len() - 1
    }
}

impl<T: Copy, A: Allocator + Copy> TwoDimArr<T, A> {
    #[inline]
    pub fn push<Arr: IntoIterator<Item = T>>(&mut self, arr: Arr) {
        self.data.extend(arr);
        self.separators.push(self.data.len());
    }
}

impl<T, A: Allocator + Copy> Index<usize> for TwoDimArr<T, A> {
    type Output = [T];

    #[inline]
    fn index(&self, idx: usize) -> &Self::Output {
        let start = self.separators[idx];
        let end = self.separators[idx + 1];
        &self.data[start..end]
    }
}

pub struct TwoDimArrIter<'a, T, A: Allocator + Copy> {
    arr: &'a TwoDimArr<T, A>,
    idx: usize,
}

impl<'a, T, A: Allocator + Copy> Iterator for TwoDimArrIter<'a, T, A> {
    type Item = &'a [T];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + 1 < self.arr.separators.len() {
            let start = self.arr.separators[self.idx];
            let end = self.arr.separators[self.idx + 1];
            self.idx += 1;
            Some(&self.arr.data[start..end])
        } else {
            None
        }
    }
}
