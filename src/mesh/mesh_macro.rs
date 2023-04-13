#[macro_export]
macro_rules! build_connect_info {
	() => {
		fn n_vertices(&self) -> usize {
			return self.n_vertices;
		}
	};
}