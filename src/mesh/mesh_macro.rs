#[macro_export]
macro_rules! build_connect_info {
    () => {
        /// get bump
        #[inline(always)]
        fn bump(&self) -> &'b Bump {
            self.v_halfedge_arr.bump()
        }
        /// the number of vertices
        #[inline(always)]
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
        fn vertex(&self, vid: VertexId) -> VertexIter<'_, 'b, Self> {
            VertexIter::new(vid, self)
        }

        /// start halfedge iterator
        #[inline(always)]
        fn halfedge(&self, hid: HalfedgeId) -> HalfedgeIter<'_, 'b, Self> {
            HalfedgeIter::new(hid, self)
        }

        /// start edge iterator
        #[inline(always)]
        fn edge(&self, eid: EdgeId) -> super::EdgeIter<'_, 'b, Self> {
            EdgeIter::new(eid, self)
        }

        /// start face iterator
        #[inline(always)]
        fn face(&self, fid: FaceId) -> FaceIter<'_, 'b, Self> {
            FaceIter::new(fid, self)
        }

        #[inline(always)]
        fn vertices(&self) -> VertexIter<'_, 'b, Self> {
            VertexIter::new(VertexId::from(0), self)
        }

        #[inline(always)]
        fn halfedges(&self) -> HalfedgeIter<'_, 'b, Self> {
            HalfedgeIter::new(HalfedgeId::from(0), self)
        }

        #[inline(always)]
        fn edges(&self) -> EdgeIter<'_, 'b, Self> {
            EdgeIter::new(EdgeId::from(0), self)
        }

        #[inline(always)]
        fn faces(&self) -> FaceIter<'_, 'b, Self> {
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
        #[inline]
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

        // two vertices of the edge
        #[inline]
        fn e_vertices(&self, eid: EdgeId) -> [VertexId; 2] {
            let hid = self.e_halfedge(eid);
            [self.he_vertex(hid), self.he_tip_vertex(hid)]
        }

        /// the first halfedge of the face
        #[inline(always)]
        fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
            HalfedgeId::from(self.f_halfedge_arr[fid])
        }

        #[inline]
        fn add_vertices_data<T: 'b + Clone>(
            &mut self,
            data: Weak<RefCell<VertexData<'b, T, Self>>>,
        ) {
            self.vertices_data.push(data);
        }

        #[inline]
        fn remove_vertices_data<T: 'b + Clone>(&mut self, remove: &VertexData<'b, T, Self>) {
            self.vertices_data.retain(|data| {
                data.upgrade()
                    .map(|data| !std::ptr::eq(data.as_ptr(), remove))
                    .unwrap_or(false)
            })
        }

        #[inline]
        fn add_halfedges_data<T: 'b + Clone>(
            &mut self,
            data: Weak<RefCell<HalfedgeData<'b, T, Self>>>,
        ) {
            self.halfedges_data.push(data);
        }

        #[inline]
        fn remove_halfedges_data<T: 'b + Clone>(&mut self, remove: &HalfedgeData<'b, T, Self>) {
            self.halfedges_data.retain(|data| {
                data.upgrade()
                    .map(|data| !std::ptr::eq(data.as_ptr(), remove))
                    .unwrap_or(false)
            })
        }

        #[inline]
        fn add_edges_data<T: 'b + Clone>(&mut self, data: Weak<RefCell<EdgeData<'b, T, Self>>>) {
            self.edges_data.push(data);
        }

        #[inline]
        fn remove_edges_data<T: 'b + Clone>(&mut self, remove: &EdgeData<'b, T, Self>) {
            self.edges_data.retain(|data| {
                data.upgrade()
                    .map(|data| !std::ptr::eq(data.as_ptr(), remove))
                    .unwrap_or(false)
            })
        }

        #[inline]
        fn add_faces_data<T: 'b + Clone>(&mut self, data: Weak<RefCell<FaceData<'b, T, Self>>>) {
            self.faces_data.push(data);
        }

        #[inline]
        fn remove_faces_data<T: 'b + Clone>(&mut self, remove: &FaceData<'b, T, Self>) {
            self.faces_data.retain(|data| {
                data.upgrade()
                    .map(|data| !std::ptr::eq(data.as_ptr(), remove))
                    .unwrap_or(false)
            })
        }
    };
}
