use std::alloc::Allocator;

use hashbrown::HashMap;
use itertools::Itertools;

use crate::mesh1::element::{Edge, EdgeIter, EdgeIterMut};

use super::element::{EdgeId, ElementId, FaceId, HalfedgeExt, HalfedgeId, VertexId};

use super::mesh::{
    BaseFace, BaseHalfedge, BaseMesh, BaseVertex, ElementContainer, HasBaseMesh, Mesh, MeshCore,
};

#[derive(Default, Clone)]
pub struct HEdge<P> {
    edge: EdgeId,
    sibling: HalfedgeId,
    incoming_next: HalfedgeId,
    property: P,
}

impl<P> HEdge<P> {
    pub fn new(edge: EdgeId, sibling: HalfedgeId, incoming_next: HalfedgeId, property: P) -> Self {
        Self {
            edge,
            sibling,
            incoming_next,
            property,
        }
    }
}

impl<P> HalfedgeExt for BaseHalfedge<HEdge<P>> {
    #[inline]
    fn edge(&self) -> EdgeId {
        self.property.edge
    }

    #[inline]
    fn sibling(&self) -> HalfedgeId {
        self.property.sibling
    }

    #[inline]
    fn incoming_next(&self) -> HalfedgeId {
        self.property.incoming_next
    }
}

#[derive(Default, Clone)]
pub struct BaseEdge<P> {
    halfedge: HalfedgeId,
    property: P,
}

impl<P> Edge for BaseEdge<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }
}

impl<P> BaseEdge<P> {
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        Self { halfedge, property }
    }
}

pub struct SurfaceMesh<VP, HP, EP, FP, A: Allocator> {
    base: BaseMesh<BaseVertex<VP>, BaseHalfedge<HEdge<HP>>, BaseFace<FP>, A>,
    edges: ElementContainer<BaseEdge<EP>, EdgeId, A>,
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

                let he = mesh.halfedge_mut(hid);
                he.vertex = vid;
                he.face = fid;

                if i == 0 {
                    mesh.new_faces(1);
                    mesh.face_mut(fid).halfedge = hid;
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
                let he = mesh.halfedge_mut(hid);
                he.property.sibling = *prev_hid;
                he.property.edge = eid;
                *prev_hid = hid;
            } else {
                // This is the first time we've ever seen this edge, create a new edge object
                let new_eid = mesh.new_edges(1);
                let he = mesh.halfedge_mut(hid);
                he.property.edge = new_eid;
                he.property.sibling = HalfedgeId::default();
                mesh.edge_mut(new_eid).halfedge = hid;
                edge_history.insert(key, hid);
            }
        }
        // Complete the sibling cycle by following backwards each edge until we reach the first sibling-less entry
        for last_hid in edge_history.into_values() {
            let mut curr_hid = last_hid;
            while mesh.halfedge(curr_hid).property.sibling.valid() {
                curr_hid = mesh.halfedge(curr_hid).property.sibling;
            }
            mesh.halfedge_mut(curr_hid).property.sibling = last_hid;
        }

        let (v_in_halfedges, v_in_separators) = mesh.vertex_cycle();
        let n_vertices = mesh.n_vertices();
        for idx in 0..n_vertices {
            let (start, end) = (v_in_separators[idx], v_in_separators[idx + 1]);
            for (&ha, &hb) in v_in_halfedges[start..end].iter().circular_tuple_windows() {
                mesh.halfedge_mut(ha).property.incoming_next = hb;
            }
        }

        mesh
    }

    #[inline]
    pub fn edge(&self, eid: EdgeId) -> &BaseEdge<EP> {
        &self.edges[eid]
    }

    #[inline]
    pub fn edge_mut(&mut self, eid: EdgeId) -> &mut BaseEdge<EP> {
        &mut self.edges[eid]
    }

    #[inline]
    pub fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        self.base_mut().new_halfedges(n)
    }

    #[inline]
    pub fn new_edges(&mut self, n: usize) -> EdgeId {
        let ret = EdgeId(self.n_edges_capacity());
        let cap = *ret + n;
        self.edges.data.resize(cap, BaseEdge::default());
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
}

impl<VP, HP, EP, FP, A: Allocator> HasBaseMesh for SurfaceMesh<VP, HP, EP, FP, A> {
    type A = A;
    type VP = VP;
    type HP = HEdge<HP>;
    type FP = FP;

    fn base(
        &self,
    ) -> &BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A> {
        &self.base
    }

    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A>
    {
        &mut self.base
    }
}

impl<VP, HP, EP, FP, A: Allocator> Mesh for SurfaceMesh<VP, HP, EP, FP, A> {
    type Edge = BaseEdge<EP>;

    #[inline]
    fn n_edges(&self) -> usize {
        self.n_edges
    }

    #[inline]
    fn n_edges_capacity(&self) -> usize {
        self.edges.len()
    }

    #[inline]
    fn edges(&self) -> impl Iterator<Item = &Self::Edge> {
        self.edges.iter()
    }

    #[inline]
    fn edges_mut(&mut self) -> impl Iterator<Item = &mut Self::Edge> {
        self.edges.iter_mut()
    }

    #[inline]
    fn edge(&self, eid: EdgeId) -> &Self::Edge {
        &self.edges[eid]
    }

    #[inline]
    fn edge_mut(&mut self, eid: EdgeId) -> &mut Self::Edge {
        &mut self.edges[eid]
    }

    #[inline]
    fn edge_iter(&self, eid: EdgeId) -> super::element::EdgeIter<Self> {
        EdgeIter::new(eid, self)
    }

    #[inline]
    fn edge_iter_mut(&mut self, eid: EdgeId) -> super::element::EdgeIterMut<Self> {
        EdgeIterMut::new(eid, self)
    }

    #[inline]
    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge(hid).property.sibling
    }

    #[inline]
    fn he_incoming_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge(hid).property.incoming_next
    }

    #[inline]
    fn he_edge(&self, hid: HalfedgeId) -> EdgeId {
        self.halfedge(hid).property.edge
    }

    #[inline]
    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId {
        self.edge(eid).halfedge
    }
}

mod tests {
    #[test]
    fn test_surface_mesh() {
        use super::SurfaceMesh;
        use crate::mesh1::mesh::MeshCore;
        use crate::mesh1::element::VertexEdgesAndVertices;
        use crate::mesh1::mesh::Mesh;

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
        let incoming_halfedges = Vec::from_iter(mesh.vertex_iter(0.into()).incoming_halfedge_ids());
        debug_assert!(incoming_halfedges.len() == 6);

        let outgoing_halfedges = Vec::from_iter(mesh.vertex_iter(0.into()).outgoing_halfedge_ids());
        debug_assert!(outgoing_halfedges.len() == 6);

        let edges = Vec::from_iter(mesh.vertex_iter(0.into()).edge_ids());
        debug_assert!(edges.len() == 6);

        let vertices = Vec::from_iter(mesh.vertex_iter(0.into()).vertex_ids());
        debug_assert!(vertices.len() == 6);

        let vertices1 = Vec::from_iter(mesh.vertex_iter_mut(0.into()).vertices().map(|v| v.id));
        debug_assert!(vertices1.len() == 6);

        let edge_halfedges = Vec::from_iter(mesh.edge_iter(0.into()).halfedges().map(|he| he.id));
        debug_assert!(edge_halfedges.len() == 2);

        let face_halfedges = Vec::from_iter(mesh.face_iter(0.into()).halfedges().map(|he| he.id));
        debug_assert!(face_halfedges.len() == 3);

        let face_halfedges1 = Vec::from_iter(mesh.face_iter(0.into()).halfedges().rev().map(|he| he.from()));
        debug_assert!(face_halfedges1.len() == 3);
    }
}
