#[macro_export]
macro_rules! build_connect_info {
    () => {
        /// the number of vertices
        fn n_vertices(&self) -> usize {
            return self.n_vertices;
        }
        /// the number of halfedges
        #[inline(always)]
        fn n_halfedges(&self) -> usize {
            return self.n_halfedges;
        }
        /// the number of faces
        #[inline(always)]
        fn n_faces(&self) -> usize {
            return self.n_faces;
        }

        /// the capacity of vertices
        #[inline(always)]
        fn n_vertices_capacity(&self) -> usize {
            return self.v_halfedge_arr.len();
        }
        /// the capacity of halfedges
        #[inline(always)]
        fn n_halfedges_capacity(&self) -> usize {
            return self.he_next_arr.len();
        }

        /// the capacity of faces
        #[inline(always)]
        fn n_faces_capacity(&self) -> usize {
            return self.f_halfedge_arr.len();
        }

        /// vertex id is valid
        #[inline(always)]
        fn vertex_is_valid(&self, vid: VertexId) -> bool {
            self.v_halfedge_arr[vid].valid()
        }

        /// halfedge id is valid
        #[inline(always)]
        fn halfedge_is_valid(&self, hid: HalfedgeId) -> bool {
            self.he_next_arr[hid].valid()
        }

        /// face id is valid
        #[inline(always)]
        fn face_is_valid(&self, fid: FaceId) -> bool {
            self.f_halfedge_arr[fid].valid()
        }

        /// start vertex iterator
        #[inline(always)]
        fn vertex<'a>(&'a self, vid: VertexId) -> VertexIter<'a, Self> {
            VertexIter::new(vid, self)
        }

        /// start halfedge iterator
        #[inline(always)]
        fn halfedge<'a>(&'a self, hid: HalfedgeId) -> HalfedgeIter<'a, Self> {
            HalfedgeIter::new(hid, self)
        }

        /// start edge iterator
        #[inline(always)]
        fn edge<'a>(&'a self, eid: EdgeId) -> super::EdgeIter<'a, Self> {
            EdgeIter::new(eid, self)
        }

        /// start face iterator
        #[inline(always)]
        fn face<'a>(&'a self, fid: FaceId) -> FaceIter<'a, Self> {
            FaceIter::new(fid, self)
        }

        #[inline(always)]
        fn vertices(&self) -> VertexIter<'_, Self> {
            VertexIter::new(VertexId::from(0), self)
        }

        #[inline(always)]
        fn halfedges(&self) -> HalfedgeIter<'_, Self> {
            HalfedgeIter::new(HalfedgeId::from(0), self)
        }

        #[inline(always)]
        fn edges(&self) -> EdgeIter<'_, Self> {
            EdgeIter::new(EdgeId::from(0), self)
        }

        #[inline(always)]
        fn faces(&self) -> FaceIter<'_, Self> {
            FaceIter::new(FaceId::from(0), self)
        }
        /// the halfedge starting from this vertex
        #[inline(always)]
        fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
            HalfedgeId::from(self.v_halfedge_arr[vid])
        }

        /// the start vertex of the halfedge
        #[inline(always)]
        fn he_vertex(&self, hid: HalfedgeId) -> VertexId {
            VertexId::from(self.he_vertex_arr[hid])
        }

        /// the end vertex of the halfedge
        #[inline(always)]
        fn he_tip_vertex(&self, hid: HalfedgeId) -> VertexId {
            VertexId::from(self.he_vertex_arr[self.he_next_arr[hid]])
        }

        /// the next halfedge of the halfedge
        #[inline(always)]
        fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
            HalfedgeId::from(self.he_next_arr[hid])
        }

        /// the previous halfedge of the halfedge
        fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
            let mut curr = hid;
            loop {
                let next = self.he_next(curr);
                if next == hid {
                    return curr;
                }
                curr = next;
            }
        }

        /// the first halfedge of the face
        #[inline(always)]
        fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
            HalfedgeId::from(self.f_halfedge_arr[fid])
        }

        #[inline]
        fn add_halfedges_data<T: Default + Clone + 'static>(
            &mut self,
            data: Weak<RefCell<HalfedgeData<T, Self>>>,
        ) {
            self.halfedges_data.push(data);
        }

        #[inline]
        fn remove_halfedges_data<T: Default + Clone>(&mut self, remove: &HalfedgeData<T, Self>) {
            self.halfedges_data.retain(|data| {
                data.upgrade()
                    .map(|data| !std::ptr::eq(data.as_ptr(), remove))
                    .unwrap_or(false)
            })
        }
    };
}
