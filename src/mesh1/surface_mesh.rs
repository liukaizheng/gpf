use std::alloc::Allocator;
use std::ptr::NonNull;

use hashbrown::HashMap;
use itertools::Itertools;

use crate::mesh1::element::{
    Edge, EdgeData, EdgeHalfedge, EdgeMut, FaceData, HalfedgeData, HalfedgeMut,
};

use super::element::{EdgeId, ElementId, FaceId, HalfedgeDataExt, HalfedgeId, VertexId};

use super::mesh::{
    BaseFaceData, BaseHalfedgeData, BaseMesh, BaseVertexData, ElementContainer, HasBaseMesh, Mesh,
    MeshCore, edge_from_vertices,
};

#[derive(Default, Clone)]
pub struct ExtendedHalfedge<P> {
    edge: EdgeId,
    sibling: HalfedgeId,
    incoming_next: HalfedgeId,
    pub property: P,
}

impl<P> ExtendedHalfedge<P> {
    pub fn new(edge: EdgeId, sibling: HalfedgeId, incoming_next: HalfedgeId, property: P) -> Self {
        Self {
            edge,
            sibling,
            incoming_next,
            property,
        }
    }
}

impl<P> HalfedgeDataExt for BaseHalfedgeData<ExtendedHalfedge<P>> {
    #[inline]
    fn edge(&self) -> EdgeId {
        self.property.edge
    }

    #[inline]
    fn set_edge(&mut self, edge: EdgeId) {
        self.property.edge = edge;
    }

    #[inline]
    fn sibling(&self) -> HalfedgeId {
        self.property.sibling
    }

    #[inline]
    fn set_sibling(&mut self, sibling: HalfedgeId) {
        self.property.sibling = sibling;
    }

    #[inline]
    fn incoming_next(&self) -> HalfedgeId {
        self.property.incoming_next
    }

    #[inline]
    fn set_incoming_next(&mut self, incoming_next: HalfedgeId) {
        self.property.incoming_next = incoming_next;
    }
}

#[derive(Default, Clone)]
pub struct BaseEdgeData<P> {
    halfedge: HalfedgeId,
    property: P,
}

impl<P> EdgeData for BaseEdgeData<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }
}

impl<P> BaseEdgeData<P> {
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        Self { halfedge, property }
    }
}

pub struct SurfaceMesh<VP, HP, EP, FP, A: Allocator> {
    base: BaseMesh<BaseVertexData<VP>, BaseHalfedgeData<ExtendedHalfedge<HP>>, BaseFaceData<FP>, A>,
    edges: ElementContainer<BaseEdgeData<EP>, EdgeId, A>,
    n_edges: usize,
}

impl<
    VP: Default + Clone,
    HP: Default + Clone,
    EP: Default + Clone,
    FP: Default + Clone,
    A: Allocator + Copy,
> SurfaceMesh<VP, HP, EP, FP, A>
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
            n_edges: 0,
        };

        for (fid, polygon) in polygons.into_iter().enumerate() {
            let fid = fid.into();
            let mut first_hid = HalfedgeId::default();
            let mut prev_hid = first_hid;
            let mut prev_vid = VertexId::default();

            for (i, &b) in polygon.as_ref().iter().enumerate() {
                let vid = VertexId::from(b);
                if vid.valid() {
                    mesh.base.v_min_reserve(vid);
                }
                let hid = mesh.new_halfedges(1);

                let he = mesh.halfedge_data_mut(hid);
                he.vertex = vid;
                he.face = fid;

                if i == 0 {
                    mesh.new_faces(1);
                    mesh.face_data_mut(fid).halfedge = hid;
                    first_hid = hid;
                } else {
                    mesh.set_v_halfedge(prev_vid, hid);
                    mesh.connect_halfedges(prev_hid, hid);
                }
                prev_vid = vid;
                prev_hid = hid;
            }
            mesh.set_v_halfedge(prev_vid, first_hid);
            mesh.connect_halfedges(prev_hid, first_hid);
        }
        mesh.base_mut().recount_n_vertices();

        let mut edge_history = HashMap::<(VertexId, VertexId), HalfedgeId>::new();
        // build edge
        for hid in 0..mesh.n_halfedges_capacity() {
            let hid = hid.into();
            let [va, vb] = mesh.he_vertices(hid);
            let key = if *va < *vb { (va, vb) } else { (vb, va) };
            if let Some(prev_hid) = edge_history.get_mut(&key) {
                // We're already seen this edge, connect to the previous halfedge incident on the edge
                let eid = mesh.he_edge(*prev_hid);
                let he = mesh.halfedge_data_mut(hid);
                he.property.sibling = *prev_hid;
                he.property.edge = eid;
                *prev_hid = hid;
            } else {
                // This is the first time we've ever seen this edge, create a new edge object
                let new_eid = mesh.new_edges(1);
                let he = mesh.halfedge_data_mut(hid);
                he.property.edge = new_eid;
                he.property.sibling = HalfedgeId::default();
                mesh.edge_data_mut(new_eid).halfedge = hid;
                edge_history.insert(key, hid);
            }
        }
        // Complete the sibling cycle by following backwards each edge until we reach the first sibling-less entry
        for last_hid in edge_history.into_values() {
            let mut curr_hid = last_hid;
            while mesh.halfedge_data(curr_hid).property.sibling.valid() {
                curr_hid = mesh.halfedge_data(curr_hid).property.sibling;
            }
            mesh.halfedge_data_mut(curr_hid).property.sibling = last_hid;
        }

        let (v_in_halfedges, v_in_separators) = mesh.vertex_cycle();
        let n_vertices = mesh.n_vertices();
        for idx in 0..n_vertices {
            let (start, end) = (v_in_separators[idx], v_in_separators[idx + 1]);
            for (&ha, &hb) in v_in_halfedges[start..end].iter().circular_tuple_windows() {
                mesh.halfedge_data_mut(ha).property.incoming_next = hb;
            }
        }

        mesh
    }

    #[inline]
    pub fn edge_data(&self, eid: EdgeId) -> &BaseEdgeData<EP> {
        &self.edges[eid]
    }

    #[inline]
    pub fn edge_data_mut(&mut self, eid: EdgeId) -> &mut BaseEdgeData<EP> {
        &mut self.edges[eid]
    }

    #[inline]
    pub fn new_vertices(&mut self, n: usize) -> VertexId {
        self.base_mut().new_vertices(n)
    }

    #[inline]
    pub fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        self.base_mut().new_halfedges(n)
    }

    #[inline]
    pub fn new_edges(&mut self, n: usize) -> EdgeId {
        let ret = EdgeId(self.n_edges_capacity());
        let cap = *ret + n;
        self.edges.data.resize(cap, BaseEdgeData::default());
        self.n_edges += n;
        ret
    }

    #[inline]
    pub fn new_faces(&mut self, n: usize) -> FaceId {
        self.base.new_faces(n)
    }

    fn vertex_cycle(&self) -> (Vec<HalfedgeId>, Vec<usize>) {
        let mut v_degree = vec![0usize; self.n_vertices_capacity()];
        for he in self.base.halfedges.iter() {
            let vid = he.vertex;
            if vid.valid() {
                v_degree[*vid] += 1;
            }
        }
        let mut vertex_separators = vec![0];
        vertex_separators.extend(v_degree.iter().scan(0, |sum, &count| {
            *sum += count;
            Some(*sum)
        }));
        let mut he_positions = vertex_separators.clone();
        let mut vertex_halfedges = vec![HalfedgeId::from(0); self.n_halfedges_capacity()];
        self.base
            .halfedges
            .iter()
            .enumerate()
            .for_each(|(hid, he)| {
                let vid = he.vertex;
                if vid.valid() {
                    let pos = he_positions[*vid];
                    vertex_halfedges[pos] = hid.into();
                    he_positions[*vid] += 1;
                }
            });
        (vertex_halfedges, vertex_separators)
    }

    pub fn split_edge(&mut self, eid: EdgeId) -> VertexId {
        let (vb, edge_hid, n_halfedges) = {
            let edge = self.edge(eid);
            let he = edge.halfedge();
            (edge.halfedge().data.vertex, he.id, edge.halfedges().count())
        };
        let start_hid = self.new_halfedges(n_halfedges);
        let new_eid = self.new_edges(1);
        let new_vid = self.new_vertices(1);
        self.vertex_data_mut(new_vid).halfedge = edge_hid;

        let mut old_prev_he: Option<HalfedgeMut<'_, Self>> = None;
        let mut new_prev_he: Option<HalfedgeMut<'_, Self>> = None;
        let mut old_first_hid = HalfedgeId::default();
        let mut new_first_hid = HalfedgeId::default();
        let mut mesh = NonNull::from_mut(self);
        unsafe {
            let mut old_edge = mesh.as_mut().edge_mut(eid);
            for (mut raw_old_he, mut raw_new_he) in old_edge
                .halfedges_mut()
                .zip(mesh.as_mut().halfedge_range_mut(start_hid, n_halfedges))
            {
                let new_hid = raw_new_he.id;
                raw_new_he.data.vertex = new_vid;
                let from_vertex = raw_old_he.from_mut();
                if from_vertex.data.halfedge == raw_old_he.id {
                    from_vertex.data.halfedge = raw_new_he.id;
                }

                let mut prev_he = raw_old_he.prev_mut();
                prev_he.connect(&mut raw_new_he);
                raw_new_he.connect(&mut raw_old_he);
                raw_new_he.data.set_face(raw_old_he.data.face);

                let [old_he, new_he] = if from_vertex.id != vb {
                    raw_new_he.data.property.edge = new_eid;
                    [raw_old_he, raw_new_he]
                } else {
                    raw_old_he.data.property.edge = new_eid;
                    raw_new_he.data.property.edge = eid;
                    [raw_new_he, raw_old_he]
                };

                if !old_first_hid.valid() {
                    old_first_hid = old_he.id;
                    new_first_hid = new_he.id;
                } else {
                    let old_prev_he = old_prev_he.as_mut().unwrap_unchecked();
                    let new_prev_he = new_prev_he.as_mut().unwrap_unchecked();
                    old_prev_he.data.property.sibling = old_he.id;
                    new_prev_he.data.property.sibling = new_he.id;
                    if old_prev_he.data.vertex == new_vid {
                        old_prev_he.data.property.incoming_next = new_hid;
                    } else {
                        debug_assert!(new_prev_he.data.vertex == new_vid);
                        new_prev_he.data.property.incoming_next = new_hid;
                    }
                }
                old_prev_he.replace(old_he);
                new_prev_he.replace(new_he);
            }
            old_edge.data.halfedge = old_first_hid;
            self.edge_data_mut(new_eid).halfedge = new_first_hid;

            let old_prev_he = old_prev_he.unwrap_unchecked();
            let new_prev_he = new_prev_he.unwrap_unchecked();
            if old_prev_he.data.vertex == new_vid {
                old_prev_he.data.property.incoming_next = start_hid;
            } else {
                debug_assert!(new_prev_he.data.vertex == new_vid);
                new_prev_he.data.property.incoming_next = start_hid;
            }

            old_prev_he.data.property.sibling = old_first_hid;
            new_prev_he.data.property.sibling = new_first_hid;
        }

        new_vid
    }

    pub fn split_face(&mut self, fid: FaceId, va: VertexId, vb: VertexId) -> HalfedgeId {
        let new_start_hid = self.new_halfedges(2);
        let new_eid = self.new_edges(1);
        let new_fid = self.new_faces(1);

        let mesh_ptr = self as *mut Self;
        let mut left_last_he = unsafe {
            (*mesh_ptr)
                .vertex_mut(va)
                .incoming_halfedges_mut()
                .find(|he| he.data.face == fid)
                .unwrap_unchecked()
        };

        let mut right_last_he = unsafe {
            (*mesh_ptr)
                .vertex_mut(vb)
                .incoming_halfedges_mut()
                .find(|he| he.data.face == fid)
                .unwrap_unchecked()
        };

        let mut left_first_he = right_last_he.next_mut();
        let mut right_first_he = left_last_he.next_mut();

        let mut first_he = unsafe { (*mesh_ptr).halfedge_mut(new_start_hid) };
        let mut second_he = unsafe { (*mesh_ptr).halfedge_mut((*new_start_hid + 1).into()) };

        right_last_he.connect(&mut first_he);
        first_he.connect(&mut right_first_he);
        left_last_he.connect(&mut second_he);
        second_he.connect(&mut left_first_he);

        first_he.data.set_vertex(va);
        first_he.data.set_edge(new_eid);
        first_he.data.set_face(fid);

        second_he.data.set_vertex(vb);
        second_he.data.set_edge(new_eid);
        second_he.data.set_face(new_fid);

        self.set_e_halfedge(new_eid, second_he.id);

        left_last_he.insert_incoming_next(&mut first_he);
        right_last_he.insert_incoming_next(&mut second_he);

        {
            let mut curr_he = left_first_he;
            loop {
                curr_he.data.face = new_fid;
                curr_he = curr_he.next_mut();
                if curr_he.id == second_he.id {
                    break;
                }
            }
        }

        self.face_mut(fid).data.set_halfedge(first_he.id);
        self.face_mut(new_fid).data.set_halfedge(second_he.id);
        second_he.id
    }

    #[inline]
    fn set_e_halfedge(&mut self, eid: EdgeId, hid: HalfedgeId) {
        self.edges[eid].halfedge = hid;
    }
}

impl<VP, HP, EP, FP, A: Allocator> HasBaseMesh for SurfaceMesh<VP, HP, EP, FP, A> {
    type A = A;
    type VP = VP;
    type HP = ExtendedHalfedge<HP>;
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

impl<VP, HP, EP, FP, A: Allocator> SurfaceMesh<VP, HP, EP, FP, A> {
    #[inline]
    pub fn edge_data_range(
        &self,
        eid: EdgeId,
        count: usize,
    ) -> impl Iterator<Item = &BaseEdgeData<EP>> {
        let end = (*eid + count).min(self.edges.len());
        self.edges.range(*eid, end)
    }

    #[inline]
    pub fn edge_data_range_mut(
        &mut self,
        eid: EdgeId,
        count: usize,
    ) -> impl Iterator<Item = &mut BaseEdgeData<EP>> {
        let end = (*eid + count).min(self.edges.len());
        self.edges.range_mut(*eid, end)
    }
}

impl<VP, HP, EP, FP, A: Allocator> Mesh for SurfaceMesh<VP, HP, EP, FP, A> {
    type EdgeData = BaseEdgeData<EP>;

    #[inline]
    fn n_edges(&self) -> usize {
        self.n_edges
    }

    #[inline]
    fn n_edges_capacity(&self) -> usize {
        self.edges.len()
    }

    #[inline]
    fn edge_datum(&self) -> impl Iterator<Item = &Self::EdgeData> {
        self.edges.iter()
    }

    #[inline]
    fn edge_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::EdgeData> {
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

    fn edge_range(&'_ self, eid: EdgeId, count: usize) -> impl Iterator<Item = Edge<'_, Self>> {
        self.edge_data_range(eid, count)
            .zip(*eid..*eid + count)
            .map(|(data, eid)| Edge::new_with_data(eid.into(), data, self))
    }

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
    fn edge(&'_ self, eid: EdgeId) -> super::element::Edge<'_, Self> {
        Edge::new(eid, self)
    }

    #[inline]
    fn edge_mut(&'_ mut self, eid: EdgeId) -> super::element::EdgeMut<'_, Self> {
        EdgeMut::new(eid, self)
    }

    #[inline]
    fn edges(&'_ self) -> impl Iterator<Item = Edge<'_, Self>> {
        self.edge_datum()
            .zip(0..self.n_edges_capacity())
            .filter_map(|(data, eid)| {
                if data.halfedge.valid() {
                    Some(Edge::new_with_data(eid.into(), data, self))
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
                if data.halfedge.valid() {
                    Some(EdgeMut::new(eid.into(), &mut *mesh_ptr))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge_data(hid).property.sibling
    }

    #[inline]
    fn he_incoming_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge_data(hid).property.incoming_next
    }

    #[inline]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        self.halfedge_data(hid).property.edge
    }

    #[inline]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        self.edge_data(eid).halfedge
    }

    #[inline]
    fn e_from_vertices(&self, va: VertexId, vb: VertexId) -> EdgeId {
        edge_from_vertices(self, va, vb)
    }
}

mod tests {

    #[test]
    fn test_surface_mesh() {
        use super::SurfaceMesh;
        use crate::mesh1::mesh::Mesh;
        use crate::mesh1::mesh::MeshCore;

        let mut mesh = SurfaceMesh::<(), (), (), (), _>::new_in(
            vec![
                vec![0, 1, 2],
                vec![0, 2, 3],
                vec![0, 3, 1],
                vec![0, 4, 5],
                vec![0, 5, 6],
                vec![0, 6, 4],
            ],
            std::alloc::Global,
        );
        debug_assert!(mesh.n_vertices() == 7);
        println!("the size of mesh is {:?}", std::mem::size_of_val(&mesh));
        let incoming_halfedges =
            Vec::from_iter(mesh.vertex(0.into()).incoming_halfedges().map(|he| he.id));
        debug_assert!(incoming_halfedges.len() == 6);

        let outgoing_halfedges =
            Vec::from_iter(mesh.vertex(0.into()).outgoing_halfedges().map(|he| he.id));
        debug_assert!(outgoing_halfedges.len() == 6);

        let edges = Vec::from_iter(mesh.vertex(0.into()).edges().map(|e| e.id));
        debug_assert!(edges.len() == 6);

        let vertices = Vec::from_iter(mesh.vertex(0.into()).vertices().map(|v| v.id));
        debug_assert!(vertices.len() == 6);

        let incoming_halfedges1 =
            Vec::from_iter(mesh.vertex(6.into()).incoming_halfedges().map(|he| he.id));
        debug_assert!(incoming_halfedges1.len() == 2);

        let outgoing_halfedges1 =
            Vec::from_iter(mesh.vertex(6.into()).outgoing_halfedges().map(|he| he.id));
        debug_assert!(outgoing_halfedges1.len() == 2);

        let edge1 = Vec::from_iter(mesh.vertex(6.into()).edges().map(|e| e.id));
        debug_assert!(edge1.len() == 3);

        let vertices1 = Vec::from_iter(mesh.vertex_mut(6.into()).vertices().map(|v| v.id));
        debug_assert!(vertices1.len() == 3);

        let edge_halfedges = Vec::from_iter(mesh.edge(0.into()).halfedges().map(|he| he.id));
        debug_assert!(edge_halfedges.len() == 2);

        let face_halfedges = Vec::from_iter(mesh.face(0.into()).halfedges().map(|he| he.id));
        debug_assert!(face_halfedges.len() == 3);

        let face_halfedges1 =
            Vec::from_iter(mesh.face(0.into()).halfedges().rev().map(|he| he.from().id));
        debug_assert!(face_halfedges1.len() == 3);
    }

    #[test]
    fn test_split_edge() {
        use super::SurfaceMesh;
        use crate::mesh1::element::VertexId;
        use crate::mesh1::element::{HalfedgeDataExt, HalfedgeNavigation};
        use crate::mesh1::mesh::Mesh;
        use crate::mesh1::mesh::MeshCore;
        use itertools::Itertools;
        use std::collections::HashSet;

        let mut mesh = SurfaceMesh::<(), (), (), (), _>::new_in(
            vec![
                vec![0, 1, 2],
                vec![0, 1, 3],
                vec![1, 0, 4],
                vec![0, 1, 5],
                vec![1, 0, 6],
            ],
            std::alloc::Global,
        );
        let eid = mesh.e_from_vertices(0.into(), 1.into());
        let edge_n_halfedges = mesh.edge(eid).halfedges().count();
        let new_vid = mesh.split_edge(eid);
        for face in mesh.faces() {
            let halfedges = face.halfedges().collect::<Vec<_>>();
            assert_eq!(halfedges.len(), 4);
            for (he1, he2) in halfedges.iter().circular_tuple_windows() {
                debug_assert!(he1.data.next == he2.id);
                debug_assert!(he2.data.prev == he1.id);
            }
            let mut vertices = HashSet::<VertexId, _>::new();
            for he in halfedges.iter() {
                let v = he.data.vertex;
                debug_assert!(he.data.face == face.id);
                vertices.insert(v);
            }
            debug_assert!(vertices.len() == 4);
        }

        let old_halfedges = mesh.edge(eid).halfedges().collect_vec();
        debug_assert!(old_halfedges.len() == edge_n_halfedges);
        for he in old_halfedges {
            debug_assert!(he.data.edge() == eid);
        }
        let new_edge = mesh.vertex(new_vid).halfedge().prev().edge();
        let new_eid = new_edge.id;
        let new_halfedges = new_edge.halfedges().collect_vec();
        debug_assert!(new_halfedges.len() == edge_n_halfedges);
        for he in new_halfedges {
            debug_assert!(he.data.edge() == new_eid);
        }
    }
}
