use std::{alloc::Allocator, collections::VecDeque};

use hashbrown::HashMap;
use itertools::Itertools;

use super::{
    element::{Edge, EdgeId, EdgeMut, ElementId, FaceId, HalfedgeId, VertexId},
    mesh::{
        BaseFaceData, BaseHalfedgeData, BaseMesh, BaseVertexData, ElementContainer, HasBaseMesh,
        Mesh, MeshCore, edge_from_vertices, edge_vertices, halfedge_from_oppo_vertex,
    },
};

pub struct ManifoldMesh<VP, HP, EP, FP, A: Allocator = std::alloc::Global> {
    base: BaseMesh<BaseVertexData<VP>, BaseHalfedgeData<HP>, BaseFaceData<FP>, A>,
    edges: ElementContainer<EP, EdgeId, A>,
    edges_cache: VecDeque<EdgeId, A>,
}

impl<VP, HP, EP, FP, A: Allocator> ManifoldMesh<VP, HP, EP, FP, A> {
    #[inline(always)]
    pub fn he_is_boundary(&self, hid: HalfedgeId) -> bool {
        !self.he_face(hid).valid()
    }

    #[inline]
    pub fn v_halfedge(&self, vid: VertexId) -> HalfedgeId {
        self.vertex_data(vid).halfedge
    }

    #[inline]
    fn he_is_valid(&self, hid: HalfedgeId) -> bool {
        self.halfedge_data(hid).vertex.valid()
    }
}

impl<
    VP: Default + Clone,
    HP: Default + Clone,
    EP: Default + Clone,
    FP: Default + Clone,
    A: Allocator + Copy,
> ManifoldMesh<VP, HP, EP, FP, A>
{
    pub fn new_in<U, T>(polygons: T, alloc: A) -> Self
    where
        U: AsRef<[usize]>,
        T: IntoIterator<Item = U>,
    {
        let base = BaseMesh::new_in(alloc);
        let mut mesh = Self {
            base,
            edges: ElementContainer::new_in(alloc),
            edges_cache: VecDeque::new_in(alloc),
        };

        let mut edge_map = HashMap::<(usize, usize), HalfedgeId>::new();

        for polygon in polygons.into_iter() {
            let fid = mesh.new_face();
            let mut first_hid = HalfedgeId::default();
            let mut prev_hid = HalfedgeId::default();

            for (&a, &b) in polygon.as_ref().iter().circular_tuple_windows() {
                let key = if a < b { (a, b) } else { (b, a) };
                let va = VertexId::from(a);
                let vb = VertexId::from(b);

                if va.valid() {
                    mesh.base.v_min_reserve(va);
                }
                if vb.valid() {
                    mesh.base.v_min_reserve(vb);
                }

                let hid = match edge_map.entry(key) {
                    hashbrown::hash_map::Entry::Occupied(entry) => {
                        let hid = mesh.he_twin(*entry.get());
                        *entry.into_mut() = hid;
                        hid
                    }
                    hashbrown::hash_map::Entry::Vacant(entry) => {
                        let hid = mesh.new_edge();
                        entry.insert(hid);
                        hid
                    }
                };

                if va.valid() {
                    mesh.set_v_halfedge(va, hid);
                }

                mesh.set_he_vertex(hid, vb);
                mesh.set_he_vertex(mesh.he_twin(hid), va);
                mesh.halfedge_data_mut(hid).face = fid;

                if !first_hid.valid() {
                    mesh.face_data_mut(fid).halfedge = hid;
                    first_hid = hid;
                } else {
                    mesh.connect_halfedges(prev_hid, hid);
                }
                prev_hid = hid;
            }

            mesh.connect_halfedges(prev_hid, first_hid);
        }

        mesh.base.recount_n_vertices();

        let mut edge_visited = vec![false; mesh.n_edges_capacity()];
        for twin_hid in edge_map.into_values() {
            // In the construction above, edges visited only once keep the even halfedge id in the
            // map value. Edges visited twice end up with the odd halfedge id (the second face
            // uses the twin).
            if twin_hid.0 & 1 != 0 {
                continue;
            }
            let eid = mesh.he_edge(twin_hid);
            if edge_visited[eid.0] {
                continue;
            }
            edge_visited[eid.0] = true;

            let first_hid = mesh.he_twin(twin_hid);
            debug_assert!(mesh.he_is_boundary(first_hid));
            let mut curr_hid = first_hid;
            loop {
                let mut prev_hid = curr_hid;
                loop {
                    prev_hid = mesh.he_next(mesh.he_twin(prev_hid));
                    let twin_prev_hid = mesh.he_twin(prev_hid);
                    if mesh.he_is_boundary(twin_prev_hid) {
                        prev_hid = twin_prev_hid;
                        break;
                    }
                }

                let vid = mesh.he_to(prev_hid);
                if vid.valid() {
                    mesh.set_v_halfedge(vid, curr_hid);
                }

                edge_visited[mesh.he_edge(prev_hid).0] = true;
                mesh.connect_halfedges(prev_hid, curr_hid);
                curr_hid = prev_hid;
                if curr_hid == first_hid {
                    break;
                }
            }
        }

        mesh
    }

    #[inline]
    pub fn new_vertices(&mut self, n: usize) -> VertexId {
        self.base.new_vertices(n)
    }

    #[inline]
    pub fn new_face(&mut self) -> FaceId {
        self.base.new_faces(1)
    }

    #[inline]
    pub fn new_edge(&mut self) -> HalfedgeId {
        debug_assert!(
            self.base.n_halfedges_capacity() % 2 == 0,
            "halfedges capacity must stay even for implicit twins"
        );

        if let Some(eid) = self.edges_cache.pop_back() {
            debug_assert!(eid.valid());
            debug_assert!(!self.e_is_valid(eid), "cached edges must be invalid");

            let hid = HalfedgeId::from(eid.0 << 1);
            debug_assert!(hid.valid());
            debug_assert!(hid.0 & 1 == 0, "edges must start at even halfedge ids");

            let twin_hid = self.he_twin(hid);
            *self.halfedge_data_mut(hid) = BaseHalfedgeData::default();
            *self.halfedge_data_mut(twin_hid) = BaseHalfedgeData::default();
            self.edges[eid] = EP::default();

            self.base.n_halfedges = self.base.n_halfedges.saturating_add(2);
            return hid;
        }

        let hid = self.base.new_halfedges(2);
        debug_assert!(hid.valid());
        debug_assert!(hid.0 & 1 == 0, "new edges must start at even halfedge ids");

        let eid = self.he_edge(hid);
        let len = eid.0 + 1;
        if self.edges.len() < len {
            self.edges.data.resize(len, EP::default());
        }

        hid
    }

    #[inline]
    pub fn remove_vertex(&mut self, vid: VertexId) {
        if !vid.valid() {
            return;
        }
        self.vertex_data_mut(vid).halfedge = HalfedgeId::default();
        self.base.n_vertices = self.base.n_vertices.saturating_sub(1);
    }

    #[inline]
    pub fn remove_edge(&mut self, eid: EdgeId) {
        if !eid.valid() || !self.e_is_valid(eid) {
            return;
        }
        let hid = self.e_halfedge(eid);
        let twin_hid = self.he_twin(hid);

        self.halfedge_data_mut(hid).vertex = VertexId::default();
        self.halfedge_data_mut(twin_hid).vertex = VertexId::default();
        self.edges[eid] = EP::default();
        self.edges_cache.push_back(eid);

        self.base.n_halfedges = self.base.n_halfedges.saturating_sub(2);
    }

    pub fn remove_face(&mut self, fid: FaceId) {
        if !fid.valid() {
            return;
        }

        let first_hid = self.f_halfedge(fid);
        let prev_first_hid = self.he_prev(first_hid);
        let mut curr_hid = first_hid;
        loop {
            self.halfedge_data_mut(curr_hid).face = FaceId::default();

            let next_hid = if curr_hid == prev_first_hid {
                first_hid
            } else {
                self.he_next(curr_hid)
            };
            let rev_next_hid = self.he_twin(curr_hid);
            let rev_curr_hid = self.he_twin(next_hid);

            match [
                self.he_is_boundary(rev_next_hid),
                self.he_is_boundary(rev_curr_hid),
            ] {
                [true, true] => {
                    let vid = self.he_to(curr_hid);
                    let vh = self.v_halfedge(vid);

                    if vh == rev_next_hid && self.he_prev(rev_next_hid) == rev_curr_hid {
                        self.remove_vertex(vid);
                    } else {
                        let prev_rev_next_hid = self.he_prev(rev_next_hid);
                        debug_assert!(self.he_is_valid(prev_rev_next_hid));
                        let next_rev_curr_hid = self.he_next(rev_curr_hid);
                        debug_assert!(self.he_is_valid(next_rev_curr_hid));
                        self.connect_halfedges(prev_rev_next_hid, next_rev_curr_hid);
                        if vid.valid() {
                            self.set_v_halfedge(vid, next_rev_curr_hid);
                        }
                    }
                }
                [true, false] => {
                    self.connect_halfedges(self.he_prev(rev_next_hid), next_hid);
                    let vid = self.he_to(curr_hid);
                    if vid.valid() {
                        self.set_v_halfedge(vid, next_hid);
                    }
                }
                [false, true] => {
                    let next_rev_curr_hid = self.he_next(rev_curr_hid);
                    self.connect_halfedges(curr_hid, next_rev_curr_hid);
                    let vid = self.he_to(curr_hid);
                    if vid.valid() {
                        self.set_v_halfedge(vid, next_rev_curr_hid);
                    }
                }
                [false, false] => {
                    let vid = self.he_to(curr_hid);
                    if vid.valid() {
                        self.set_v_halfedge(vid, next_hid);
                    }
                }
            }

            if self.he_is_boundary(rev_next_hid) {
                self.remove_edge(self.he_edge(curr_hid));
            }

            curr_hid = next_hid;
            if curr_hid == first_hid {
                break;
            }
        }

        self.face_data_mut(fid).halfedge = HalfedgeId::default();
        self.base.n_faces = self.base.n_faces.saturating_sub(1);
    }

    /// Adds a new face that is bounded by the given halfedges.
    ///
    /// The input halfedges are assumed to be oriented consistently along the new face boundary.
    pub fn new_face_by_halfedges(&mut self, halfedges: &[HalfedgeId]) -> FaceId {
        for (&ha_twin, &hb_twin) in halfedges.iter().rev().circular_tuple_windows() {
            let ha = self.he_twin(ha_twin);
            let hb = self.he_twin(hb_twin);

            match [self.he_is_boundary(ha), self.he_is_boundary(hb)] {
                [true, true] => {
                    let vid = self.he_to(ha);
                    let vh = if vid.valid() {
                        self.v_halfedge(vid)
                    } else {
                        HalfedgeId::default()
                    };

                    if vid.valid() {
                        self.set_v_halfedge(vid, hb);
                    }
                    if vh.valid() {
                        let vh_prev = self.he_prev(vh);
                        if self.he_is_boundary(vh_prev) && self.he_is_boundary(vh) {
                            self.connect_halfedges(vh_prev, hb);
                            self.connect_halfedges(ha, vh);
                            continue;
                        }
                    }
                    self.connect_halfedges(ha, hb);
                }
                [true, false] => {
                    let ha_next = self.he_next(hb_twin);
                    self.connect_halfedges(ha, ha_next);
                }
                [false, true] => {
                    let hb_prev = self.he_prev(ha_twin);
                    self.connect_halfedges(hb_prev, hb);
                    let vid = self.he_to(ha);
                    if vid.valid() {
                        self.set_v_halfedge(vid, hb);
                    }
                }
                [false, false] => {}
            }
        }

        let new_fid = self.new_face();
        for (&ha, &hb) in halfedges.iter().circular_tuple_windows() {
            self.connect_halfedges(ha, hb);
            self.halfedge_data_mut(ha).face = new_fid;
        }
        self.face_data_mut(new_fid).halfedge = halfedges[0];

        new_fid
    }

    #[inline]
    pub fn he_replace(&mut self, old_hid: HalfedgeId, new_hid: HalfedgeId) {
        let va = self.he_from(old_hid);
        let prev_hid = self.he_prev(old_hid);
        let next_hid = self.he_next(old_hid);

        self.connect_halfedges(prev_hid, new_hid);
        self.connect_halfedges(new_hid, next_hid);

        let fid = self.he_face(old_hid);
        self.halfedge_data_mut(new_hid).face = fid;

        if va.valid() && self.v_halfedge(va) == old_hid {
            self.set_v_halfedge(va, new_hid);
        }

        if fid.valid() {
            self.set_f_halfedge(fid, new_hid);
            self.halfedge_data_mut(old_hid).face = FaceId::default();
        }

        if self.he_is_boundary(self.he_twin(old_hid)) {
            self.remove_edge(self.he_edge(old_hid));
        }
    }

    pub fn flip(&mut self, hid: HalfedgeId) {
        let bl_hid = self.he_next(hid);
        let br_hid = self.he_next(bl_hid);

        let twin_hid = self.he_twin(hid);
        let tr_hid = self.he_next(twin_hid);
        let tl_hid = self.he_next(tr_hid);

        let fid = self.he_face(hid);
        let twin_fid = self.he_face(twin_hid);

        let bottom_vid = self.he_to(bl_hid);
        let top_vid = self.he_to(tr_hid);
        let left_vid = self.he_to(hid);
        let right_vid = self.he_to(twin_hid);

        self.halfedge_data_mut(tr_hid).face = fid;
        self.halfedge_data_mut(bl_hid).face = twin_fid;

        self.set_he_vertex(hid, bottom_vid);
        self.set_he_vertex(twin_hid, top_vid);

        self.connect_halfedges(hid, br_hid);
        self.connect_halfedges(br_hid, tr_hid);
        self.connect_halfedges(tr_hid, hid);

        self.connect_halfedges(twin_hid, tl_hid);
        self.connect_halfedges(tl_hid, bl_hid);
        self.connect_halfedges(bl_hid, twin_hid);

        self.set_f_halfedge(fid, hid);
        self.set_f_halfedge(twin_fid, twin_hid);

        if self.v_halfedge(left_vid) == twin_hid {
            self.set_v_halfedge(left_vid, bl_hid);
        }
        if self.v_halfedge(right_vid) == hid {
            self.set_v_halfedge(right_vid, tr_hid);
        }
    }

    pub fn split_edge(&mut self, eid: EdgeId) -> VertexId {
        let hid = self.e_halfedge(eid);
        let twin_hid = self.he_twin(hid);
        let vb = self.he_to(hid);

        let new_v = self.new_vertices(1);
        let new_hid = self.new_edge();
        let new_twin_hid = self.he_twin(new_hid);

        if vb.valid() && self.v_halfedge(vb) == twin_hid {
            self.set_v_halfedge(vb, new_twin_hid);
        }
        self.set_v_halfedge(new_v, new_hid);

        let fid = self.he_face(hid);
        let twin_fid = self.he_face(twin_hid);

        self.halfedge_data_mut(new_hid).face = fid;
        self.halfedge_data_mut(new_twin_hid).face = twin_fid;

        self.set_he_vertex(hid, new_v);
        self.set_he_vertex(new_hid, vb);
        self.set_he_vertex(new_twin_hid, new_v);

        let prev_twin_hid = self.he_prev(twin_hid);
        let next_hid = self.he_next(hid);

        self.connect_halfedges(prev_twin_hid, new_twin_hid);
        self.connect_halfedges(new_twin_hid, twin_hid);

        self.connect_halfedges(hid, new_hid);
        self.connect_halfedges(new_hid, next_hid);

        new_v
    }

    pub fn split_face(&mut self, fid: FaceId, va: VertexId, vb: VertexId) -> HalfedgeId {
        let mut left_last_hid = HalfedgeId::default();
        let mut right_last_hid = HalfedgeId::default();
        for he in self.face(fid).halfedges() {
            let v = he.data.vertex;
            if v == va {
                left_last_hid = he.id;
            } else if v == vb {
                right_last_hid = he.id;
            }
        }

        debug_assert!(left_last_hid.valid());
        debug_assert!(right_last_hid.valid());

        let left_first_hid = self.he_next(right_last_hid);
        let right_first_hid = self.he_next(left_last_hid);

        let first_he = self.new_edge();
        let second_he = self.he_twin(first_he);

        self.halfedge_data_mut(first_he).vertex = va;
        self.halfedge_data_mut(second_he).vertex = vb;

        self.connect_halfedges(right_last_hid, first_he);
        self.connect_halfedges(first_he, right_first_hid);
        self.connect_halfedges(left_last_hid, second_he);
        self.connect_halfedges(second_he, left_first_hid);

        let new_f = self.new_face();
        self.halfedge_data_mut(first_he).face = fid;
        {
            let mut hid = second_he;
            loop {
                self.halfedge_data_mut(hid).face = new_f;
                hid = self.he_next(hid);
                if hid == second_he {
                    break;
                }
            }
        }

        self.face_data_mut(fid).halfedge = first_he;
        self.face_data_mut(new_f).halfedge = second_he;

        second_he
    }
}

impl<VP, HP, EP, FP, A: Allocator> HasBaseMesh for ManifoldMesh<VP, HP, EP, FP, A> {
    type A = A;
    type VP = VP;
    type HP = HP;
    type FP = FP;

    fn base(
        &self,
    ) -> &BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    > {
        &self.base
    }

    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    > {
        &mut self.base
    }
}

impl<VP, HP, EP, FP, A: Allocator> ManifoldMesh<VP, HP, EP, FP, A> {
    #[inline]
    pub fn edge_data_range(&self, eid: EdgeId, count: usize) -> impl Iterator<Item = &EP> {
        let end = (*eid + count).min(self.edges.len());
        self.edges.range(*eid, end)
    }

    #[inline]
    pub fn edge_data_range_mut(
        &mut self,
        eid: EdgeId,
        count: usize,
    ) -> impl Iterator<Item = &mut EP> {
        let end = (*eid + count).min(self.edges.len());
        self.edges.range_mut(*eid, end)
    }

    #[inline]
    fn e_is_valid(&self, eid: EdgeId) -> bool {
        let idx = eid.0 << 1;
        self.he_is_valid(idx.into()) || self.he_is_valid((idx + 1).into())
    }
}

impl<VP, HP, EP, FP, A: Allocator> Mesh for ManifoldMesh<VP, HP, EP, FP, A> {
    type EdgeData = EP;

    #[inline]
    fn n_edges(&self) -> usize {
        self.n_halfedges() >> 1
    }

    #[inline]
    fn n_edges_capacity(&self) -> usize {
        self.edges.len()
    }

    #[inline]
    fn edge_reserve(&mut self, additional: usize) {
        self.edges.reserve(additional);
        self.halfedge_reserve(additional << 1);
    }

    #[inline]
    fn edge_datum(&'_ self) -> impl Iterator<Item = &Self::EdgeData> {
        self.edges.iter()
    }

    #[inline]
    fn edge_datum_mut(&'_ mut self) -> impl Iterator<Item = &mut Self::EdgeData> {
        self.edges.iter_mut()
    }

    #[inline]
    fn edge_data(&self, eid: EdgeId) -> &Self::EdgeData {
        &self.edges[eid]
    }

    #[inline]
    fn edge_data_mut(&mut self, eid: EdgeId) -> &mut Self::EdgeData {
        &mut self.edges[eid]
    }

    #[inline]
    fn edge_range(&'_ self, eid: EdgeId, count: usize) -> impl Iterator<Item = Edge<'_, Self>> {
        self.edge_data_range(eid, count)
            .zip(*eid..*eid + count)
            .map(|(data, eid)| Edge::new_with_data(eid.into(), data, self))
    }

    #[inline]
    fn edge_range_mut(
        &'_ mut self,
        eid: EdgeId,
        count: usize,
    ) -> impl Iterator<Item = EdgeMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        self.edge_data_range_mut(eid, count)
            .zip(*eid..*eid + count)
            .map(move |(data, eid)| unsafe {
                EdgeMut::new_with_data(eid.into(), data, &mut *mesh_ptr)
            })
    }

    #[inline]
    fn edge(&'_ self, eid: EdgeId) -> Edge<'_, Self> {
        Edge::new(eid, self)
    }

    #[inline]
    fn edge_mut(&'_ mut self, eid: EdgeId) -> EdgeMut<'_, Self> {
        EdgeMut::new(eid, self)
    }

    #[inline]
    fn edges(&'_ self) -> impl Iterator<Item = Edge<'_, Self>> {
        self.edge_datum()
            .zip(0..self.n_edges_capacity())
            .filter_map(|(data, eid)| {
                let eid = EdgeId::from(eid);
                if self.e_is_valid(eid) {
                    Some(Edge::new_with_data(eid, data, self))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn edges_mut(&'_ mut self) -> impl Iterator<Item = EdgeMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        let n_edges_capacity = self.n_edges_capacity();
        self.edge_datum_mut()
            .zip(0..n_edges_capacity)
            .filter_map(move |(data, eid)| unsafe {
                let eid = EdgeId::from(eid);
                if (*mesh_ptr).e_is_valid(eid) {
                    Some(EdgeMut::new_with_data(eid, data, &mut *mesh_ptr))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        (hid.0 >> 1).into()
    }

    #[inline]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_twin(hid)
    }

    #[inline(always)]
    fn he_twin(&self, hid: HalfedgeId) -> HalfedgeId {
        (hid.0 ^ 1).into()
    }

    #[inline]
    fn he_incoming_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.he_prev(self.he_twin(hid))
    }

    #[inline]
    fn he_from_oppo_vertex(&self, fid: FaceId, vid: VertexId) -> HalfedgeId {
        halfedge_from_oppo_vertex(self, fid, vid)
    }

    #[inline]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        (eid.0 << 1).into()
    }

    #[inline]
    fn e_from_vertices(&self, va: VertexId, vb: VertexId) -> EdgeId {
        edge_from_vertices(self, va, vb)
    }

    #[inline]
    fn e_vertices(&self, eid: EdgeId) -> [VertexId; 2] {
        edge_vertices(self, eid)
    }

    #[inline]
    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId {
        self.face_data(fid).halfedge
    }
}

#[cfg(test)]
mod tests {
    use super::ManifoldMesh;
    use crate::mesh1::element::ElementId;
    use crate::mesh1::mesh::{Mesh, MeshCore};
    use crate::mesh1::{FaceId, VertexId};

    #[test]
    fn test_single_triangle_boundary_loop() {
        let mesh =
            ManifoldMesh::<(), (), (), (), _>::new_in(vec![vec![0, 1, 2]], std::alloc::Global);

        assert_eq!(mesh.n_vertices(), 3);
        assert_eq!(mesh.n_faces(), 1);
        assert_eq!(mesh.n_edges(), 3);
        assert_eq!(mesh.n_halfedges(), 6);

        assert_eq!(mesh.face(FaceId::from(0)).halfedges().count(), 3);

        let boundary_halfedges: Vec<_> = mesh
            .halfedges()
            .filter(|he| !mesh.he_face(he.id).valid())
            .map(|he| he.id)
            .collect();
        assert_eq!(boundary_halfedges.len(), 3);

        let start = boundary_halfedges[0];
        let mut curr = start;
        let mut count = 0usize;
        loop {
            count += 1;
            curr = mesh.he_next(curr);
            if curr == start {
                break;
            }
            assert!(count < 10);
        }
        assert_eq!(count, 3);

        for vid in [0usize, 1, 2].map(VertexId::from) {
            assert_eq!(mesh.vertex(vid).incoming_halfedges().count(), 2);
            assert_eq!(mesh.vertex(vid).outgoing_halfedges().count(), 2);
        }
        let halfedges = Vec::from_iter(mesh.edge(0.into()).halfedges().map(|he| he.id));
        assert_eq!(halfedges.len(), 2);
    }

    #[test]
    fn test_tetrahedron_closed() {
        let mesh = ManifoldMesh::<(), (), (), (), _>::new_in(
            vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 1], vec![1, 3, 2]],
            std::alloc::Global,
        );

        assert_eq!(mesh.n_vertices(), 4);
        assert_eq!(mesh.n_faces(), 4);
        assert_eq!(mesh.n_edges(), 6);
        assert_eq!(mesh.n_halfedges(), 12);

        assert_eq!(
            mesh.halfedges()
                .filter(|he| !mesh.he_face(he.id).valid())
                .count(),
            0
        );
        for eid in 0..mesh.n_edges_capacity() {
            assert_eq!(mesh.edge(eid.into()).halfedges().count(), 2);
        }
    }

    #[test]
    fn test_edge_reuse_cache() {
        let mut mesh = ManifoldMesh::<(), (), (), (), _>::new_in(
            Vec::<Vec<usize>>::new(),
            std::alloc::Global,
        );

        let hid0 = mesh.new_edge();
        let eid0 = mesh.he_edge(hid0);

        mesh.set_he_vertex(hid0, VertexId::from(0));
        mesh.set_he_vertex(mesh.he_twin(hid0), VertexId::from(1));

        let halfedges_capacity = mesh.n_halfedges_capacity();
        mesh.remove_edge(eid0);

        let hid1 = mesh.new_edge();
        assert_eq!(hid1, hid0);
        assert_eq!(mesh.n_halfedges_capacity(), halfedges_capacity);
    }
}
