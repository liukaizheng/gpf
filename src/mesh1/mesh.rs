use std::{
    alloc::Allocator,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use super::element::{
    EdgeId, ElementId, Face, FaceId, FaceIter, FaceIterMut, Halfedge, HalfedgeId, HalfedgeIter,
    HalfedgeIterMut, Vertex, VertexId, VertexIter, VertexIterMut,
};

pub struct ElementContainer<T, E: ElementId, A: Allocator> {
    pub(crate) data: Vec<T, A>,
    pub(crate) marker: PhantomData<E>,
}

impl<T, E: ElementId, A: Allocator> ElementContainer<T, E, A> {
    #[inline]
    pub fn new_in(alloc: A) -> Self {
        Self {
            data: Vec::new_in(alloc),
            marker: PhantomData,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    #[inline]
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }
}

impl<T, E: ElementId, A: Allocator> Index<E> for ElementContainer<T, E, A> {
    type Output = T;

    #[inline]
    fn index(&self, index: E) -> &Self::Output {
        unsafe { self.data.get_unchecked(index.index()) }
    }
}

impl<T, E: ElementId, A: Allocator> IndexMut<E> for ElementContainer<T, E, A> {
    #[inline]
    fn index_mut(&mut self, index: E) -> &mut Self::Output {
        unsafe { self.data.get_unchecked_mut(index.index()) }
    }
}

pub struct BaseMesh<VP, HP, FP, A: Allocator = std::alloc::Global> {
    pub(crate) vertices: ElementContainer<VP, VertexId, A>,
    pub(crate) halfedges: ElementContainer<HP, HalfedgeId, A>,
    pub(crate) faces: ElementContainer<FP, FaceId, A>,
    pub(crate) n_vertices: usize,
    pub(crate) n_halfedges: usize,
    pub(crate) n_faces: usize,
}

impl<VP, HP, FP, A: Allocator + Copy> BaseMesh<VP, HP, FP, A> {
    pub fn new_in(alloc: A) -> Self {
        Self {
            vertices: ElementContainer::new_in(alloc),
            halfedges: ElementContainer::new_in(alloc),
            faces: ElementContainer::new_in(alloc),
            n_vertices: 0,
            n_halfedges: 0,
            n_faces: 0,
        }
    }
}
impl<VP, HP, FP, A: Allocator> BaseMesh<VP, HP, FP, A> {
    #[inline]
    pub fn vertex(&self, vid: VertexId) -> &VP {
        &self.vertices[vid]
    }

    #[inline]
    pub fn vertex_mut(&mut self, vid: VertexId) -> &mut VP {
        &mut self.vertices[vid]
    }

    #[inline]
    pub fn halfedge(&self, hid: HalfedgeId) -> &HP {
        &self.halfedges[hid]
    }

    #[inline]
    pub fn halfedge_mut(&mut self, hid: HalfedgeId) -> &mut HP {
        &mut self.halfedges[hid]
    }

    #[inline]
    pub fn face(&self, fid: FaceId) -> &FP {
        &self.faces[fid]
    }

    #[inline]
    pub fn face_mut(&mut self, fid: FaceId) -> &mut FP {
        &mut self.faces[fid]
    }
}

impl<VP, HP, FP, A: Allocator> BaseMesh<BaseVertex<VP>, HP, FP, A> {
    #[inline]
    pub fn recount_n_vertices(&mut self) {
        self.n_vertices = self.vertices.iter().filter(|v| v.valid()).count();
    }
}

impl<VP: Default + Clone, HP, FP, A: Allocator> BaseMesh<VP, HP, FP, A> {
    #[inline]
    pub fn v_min_reserve(&mut self, vid: VertexId) {
        let len = *vid + 1;
        if self.vertices.len() < len {
            self.vertices.data.resize(len, VP::default());
        }
    }
}

impl<VP, HP: Default + Clone, FP, A: Allocator> BaseMesh<VP, HP, FP, A> {
    #[inline]
    pub fn new_halfedges(&mut self, n: usize) -> HalfedgeId {
        let ret = HalfedgeId(self.halfedges.len().into());
        self.halfedges.data.resize(*ret + n, HP::default());
        self.n_halfedges += n;
        ret
    }
}

impl<VP, HP, FP: Default + Clone, A: Allocator> BaseMesh<VP, HP, FP, A> {
    #[inline]
    pub fn new_faces(&mut self, n: usize) -> FaceId {
        let ret = FaceId(self.faces.len().into());
        self.faces.data.resize(*ret + n, FP::default());
        self.n_faces += n;
        ret
    }
}

#[derive(Default, Clone)]
pub struct BaseVertex<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> Vertex for BaseVertex<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }

    #[inline]
    fn set_halfedge(&mut self, hid: HalfedgeId) {
        self.halfedge = hid;
    }
}

impl<P> BaseVertex<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        BaseVertex { halfedge, property }
    }

    #[inline]
    fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

#[derive(Default, Clone)]
pub struct BaseHalfedge<P> {
    pub vertex: VertexId,
    pub next: HalfedgeId,
    pub prev: HalfedgeId,
    pub face: FaceId,
    pub property: P,
}

impl<P> BaseHalfedge<P> {
    #[inline]
    pub fn new(
        vertex: VertexId,
        next: HalfedgeId,
        prev: HalfedgeId,
        face: FaceId,
        property: P,
    ) -> Self {
        BaseHalfedge {
            vertex,
            next,
            prev,
            face,
            property,
        }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.vertex.valid()
    }
}

impl<P> Halfedge for BaseHalfedge<P> {
    #[inline]
    fn vertex(&self) -> VertexId {
        self.vertex
    }

    #[inline]
    fn set_vertex(&mut self, vertex: VertexId) {
        self.vertex = vertex;
    }

    #[inline]
    fn prev(&self) -> HalfedgeId {
        self.prev
    }

    #[inline]
    fn set_prev(&mut self, prev: HalfedgeId) {
        self.prev = prev;
    }

    #[inline]
    fn next(&self) -> HalfedgeId {
        self.next
    }

    #[inline]
    fn set_next(&mut self, next: HalfedgeId) {
        self.next = next;
    }

    #[inline]
    fn face(&self) -> FaceId {
        self.face
    }

    #[inline]
    fn set_face(&mut self, face: FaceId) {
        self.face = face;
    }
}

#[derive(Default, Clone)]
pub struct BaseFace<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> BaseFace<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        BaseFace { halfedge, property }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

impl<P> Face for BaseFace<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }

    #[inline]
    fn set_halfedge(&mut self, halfedge: HalfedgeId) {
        self.halfedge = halfedge;
    }
}

impl<V: Vertex, H: Halfedge, F: Face, A: Allocator> BaseMesh<V, H, F, A> {
    #[inline]
    fn he_from(&self, hid: HalfedgeId) -> VertexId {
        self.halfedge(self.halfedge(hid).prev()).vertex()
    }

    #[inline]
    fn he_to(&self, hid: HalfedgeId) -> VertexId {
        self.halfedge(hid).vertex()
    }

    #[inline]
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge(hid).prev()
    }

    #[inline]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.halfedge(hid).next()
    }

    #[inline]
    pub fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2] {
        let h2 = self.halfedge(hid);
        let h1 = self.halfedge(h2.prev());
        [h1.vertex(), h2.vertex()]
    }

    #[inline]
    fn set_v_halfedge(&mut self, v: VertexId, hid: HalfedgeId) {
        self.vertices[v].set_halfedge(hid);
    }

    #[inline]
    fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId) {
        self.faces[fid].set_halfedge(hid);
    }

    #[inline]
    fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId) {
        self.halfedges[hid].set_vertex(vid);
    }

    #[inline]
    pub fn connect_halfedges(&mut self, hid1: HalfedgeId, hid2: HalfedgeId) {
        self.halfedges[hid1].set_next(hid2);
        self.halfedges[hid2].set_prev(hid1);
    }
}

pub trait MeshCore: Sized {
    type Vertex;
    type Halfedge;
    type Face;

    fn n_vertices(&self) -> usize;
    fn n_halfedges(&self) -> usize;
    fn n_faces(&self) -> usize;

    fn n_vertices_capacity(&self) -> usize;
    fn n_halfedges_capacity(&self) -> usize;
    fn n_faces_capacity(&self) -> usize;

    fn vertex(&self, vid: VertexId) -> &Self::Vertex;
    fn halfedge(&self, hid: HalfedgeId) -> &Self::Halfedge;
    fn face(&self, fid: FaceId) -> &Self::Face;

    fn vertices(&self) -> impl Iterator<Item = &Self::Vertex>;
    fn halfedges(&self) -> impl Iterator<Item = &Self::Halfedge>;
    fn faces(&self) -> impl Iterator<Item = &Self::Face>;

    fn vertices_mut(&mut self) -> impl Iterator<Item = &mut Self::Vertex>;
    fn halfedges_mut(&mut self) -> impl Iterator<Item = &mut Self::Halfedge>;
    fn faces_mut(&mut self) -> impl Iterator<Item = &mut Self::Face>;

    fn vertex_mut(&mut self, vid: VertexId) -> &mut Self::Vertex;
    fn halfedge_mut(&mut self, hid: HalfedgeId) -> &mut Self::Halfedge;
    fn face_mut(&mut self, fid: FaceId) -> &mut Self::Face;

    fn vertex_iter(&self, vid: VertexId) -> VertexIter<Self>;
    fn halfedge_iter(&self, hid: HalfedgeId) -> HalfedgeIter<Self>;
    fn face_iter(&self, fid: FaceId) -> FaceIter<Self>;

    fn vertex_iter_mut(&mut self, vid: VertexId) -> VertexIterMut<Self>;
    fn halfedge_iter_mut(&mut self, hid: HalfedgeId) -> HalfedgeIterMut<Self>;
    fn face_iter_mut(&mut self, fid: FaceId) -> FaceIterMut<Self>;

    fn he_from(&self, hid: HalfedgeId) -> VertexId;
    fn he_to(&self, hid: HalfedgeId) -> VertexId;
    fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2];
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId;

    fn connect_halfedges(&mut self, hid1: HalfedgeId, hid2: HalfedgeId);

    fn set_v_halfedge(&mut self, v: VertexId, hid: HalfedgeId);
    fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId);
    fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId);
}

pub trait HasBaseMesh {
    type A: Allocator;
    type VP;
    type HP;
    type FP;
    fn base(
        &self,
    ) -> &BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A>;
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A>;
}

impl<VP, HP, FP, A: Allocator> HasBaseMesh
    for BaseMesh<BaseVertex<VP>, BaseHalfedge<HP>, BaseFace<FP>, A>
{
    type A = A;
    type VP = VP;
    type HP = HP;
    type FP = FP;

    #[inline]
    fn base(
        &self,
    ) -> &BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A> {
        self
    }

    #[inline]
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<BaseVertex<Self::VP>, BaseHalfedge<Self::HP>, BaseFace<Self::FP>, Self::A>
    {
        self
    }
}

impl<T: HasBaseMesh> MeshCore for T {
    type Vertex = BaseVertex<T::VP>;
    type Halfedge = BaseHalfedge<T::HP>;
    type Face = BaseFace<T::FP>;

    #[inline]
    fn n_vertices(&self) -> usize {
        self.base().n_vertices
    }
    #[inline]
    fn n_vertices_capacity(&self) -> usize {
        self.base().vertices.data.len()
    }

    #[inline]
    fn n_halfedges(&self) -> usize {
        self.base().n_halfedges
    }

    #[inline]
    fn n_faces(&self) -> usize {
        self.base().n_faces
    }

    #[inline]
    fn n_halfedges_capacity(&self) -> usize {
        self.base().halfedges.data.len()
    }

    #[inline]
    fn n_faces_capacity(&self) -> usize {
        self.base().faces.data.len()
    }

    #[inline]
    fn vertex(&self, vid: VertexId) -> &Self::Vertex {
        self.base().vertex(vid)
    }

    #[inline]
    fn halfedge(&self, hid: HalfedgeId) -> &Self::Halfedge {
        self.base().halfedge(hid)
    }

    #[inline]
    fn face(&self, fid: FaceId) -> &Self::Face {
        self.base().face(fid)
    }

    #[inline]
    fn vertex_mut(&mut self, vid: VertexId) -> &mut Self::Vertex {
        self.base_mut().vertex_mut(vid)
    }

    #[inline]
    fn halfedge_mut(&mut self, hid: HalfedgeId) -> &mut Self::Halfedge {
        self.base_mut().halfedge_mut(hid)
    }

    #[inline]
    fn face_mut(&mut self, fid: FaceId) -> &mut Self::Face {
        self.base_mut().face_mut(fid)
    }

    #[inline]
    fn vertices(&self) -> impl Iterator<Item = &Self::Vertex> {
        self.base().vertices.iter()
    }

    #[inline]
    fn halfedges(&self) -> impl Iterator<Item = &Self::Halfedge> {
        self.base().halfedges.iter()
    }

    #[inline]
    fn faces(&self) -> impl Iterator<Item = &Self::Face> {
        self.base().faces.iter()
    }

    #[inline]
    fn vertices_mut(&mut self) -> impl Iterator<Item = &mut Self::Vertex> {
        self.base_mut().vertices.iter_mut()
    }

    #[inline]
    fn halfedges_mut(&mut self) -> impl Iterator<Item = &mut Self::Halfedge> {
        self.base_mut().halfedges.iter_mut()
    }

    #[inline]
    fn faces_mut(&mut self) -> impl Iterator<Item = &mut Self::Face> {
        self.base_mut().faces.iter_mut()
    }

    #[inline]
    fn vertex_iter(&self, vid: VertexId) -> VertexIter<Self> {
        VertexIter::new(vid, self)
    }

    #[inline]
    fn halfedge_iter(&self, hid: HalfedgeId) -> HalfedgeIter<Self> {
        HalfedgeIter::new(hid, self)
    }

    #[inline]
    fn face_iter(&self, fid: FaceId) -> FaceIter<Self> {
        FaceIter::new(fid, self)
    }

    #[inline]
    fn vertex_iter_mut(&mut self, vid: VertexId) -> VertexIterMut<Self> {
        VertexIterMut::new(vid, self)
    }

    #[inline]
    fn halfedge_iter_mut(&mut self, hid: HalfedgeId) -> HalfedgeIterMut<Self> {
        HalfedgeIterMut::new(hid, self)
    }

    #[inline]
    fn face_iter_mut(&mut self, fid: FaceId) -> FaceIterMut<Self> {
        FaceIterMut::new(fid, self)
    }

    #[inline]
    fn he_from(&self, hid: HalfedgeId) -> VertexId {
        self.base().he_from(hid)
    }

    #[inline]
    fn he_to(&self, hid: HalfedgeId) -> VertexId {
        self.base().he_to(hid)
    }

    #[inline]
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId {
        self.base().he_prev(hid)
    }

    #[inline]
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId {
        self.base().he_next(hid)
    }

    #[inline]
    fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2] {
        self.base().he_vertices(hid)
    }

    #[inline]
    fn connect_halfedges(&mut self, hid1: HalfedgeId, hid2: HalfedgeId) {
        self.base_mut().connect_halfedges(hid1, hid2);
    }

    #[inline]
    fn set_v_halfedge(&mut self, vid: VertexId, hid: HalfedgeId) {
        self.base_mut().set_v_halfedge(vid, hid);
    }

    #[inline]
    fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId) {
        self.base_mut().set_f_halfedge(fid, hid);
    }

    #[inline]
    fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId) {
        self.base_mut().set_he_vertex(hid, vid);
    }
}

pub trait Mesh: MeshCore {
    type Edge;

    fn n_edges(&self) -> usize;
    fn n_edges_capacity(&self) -> usize;

    fn edges(&self) -> impl Iterator<Item = &Self::Edge>;
    fn edges_mut(&mut self) -> impl Iterator<Item = &mut Self::Edge>;

    fn edge(&self, eid: EdgeId) -> &Self::Edge;
    fn edge_mut(&mut self, eid: EdgeId) -> &mut Self::Edge;

    fn he_edge(&self, hid: HalfedgeId) -> EdgeId;

    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_incoming_next(&self, hid: HalfedgeId) -> HalfedgeId;

    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId;
}
