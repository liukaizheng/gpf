use std::alloc::Allocator;

struct TwoDimArr<T: Clone, A: Allocator + Copy = std::alloc::Global> {
    data: Vec<T, A>,
	separators: Vec<usize, A>,
}

impl <T: Clone, A: Allocator + Copy> TwoDimArr<T, A> {

	#[inline]
	pub fn new(separators: Vec<usize, A>, data: Vec<T, A>) -> Self {
		Self { data, separators }
	}

	#[inline]
	pub fn iter(&self) -> TwoDimArrIter<T, A> {
		TwoDimArrIter { arr: self, idx: 0 }
	}
}

struct TwoDimArrIter<'a, T: Clone, A: Allocator + Copy> {
	arr: &'a TwoDimArr<T, A>,
	idx: usize,
}

impl <'a, T: Clone, A: Allocator + Copy> Iterator for TwoDimArrIter<'a, T, A> {
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