use std::{
    alloc::Allocator,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use super::element::{EdgeId, ElementId, FaceId, HalfedgeId, VertexId};

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
pub struct Vertex<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> Vertex<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        Vertex { halfedge, property }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

#[derive(Default, Clone)]
pub struct Halfedge<P> {
    pub vertex: VertexId,
    pub next: HalfedgeId,
    pub prev: HalfedgeId,
    pub face: FaceId,
    pub property: P,
}

impl<P> Halfedge<P> {
    #[inline]
    pub fn new(
        vertex: VertexId,
        next: HalfedgeId,
        prev: HalfedgeId,
        face: FaceId,
        property: P,
    ) -> Self {
        Halfedge {
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

#[derive(Default, Clone)]
pub struct Face<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> Face<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        Face { halfedge, property }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

impl<VP, HP, FP, A: Allocator> BaseMesh<Vertex<VP>, Halfedge<HP>, Face<FP>, A> {
    #[inline]
    pub fn recount_n_vertices(&mut self) {
        self.n_vertices = self.vertices.iter().filter(|v| v.valid()).count();
    }


    pub fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2] {
        let h2 = self.halfedge(hid);
        let h1 = self.halfedge(h2.prev);
        [h1.vertex, h2.vertex]
    }

    #[inline]
    fn set_v_halfedge(&mut self, v: VertexId, hid: HalfedgeId) {
        self.vertices[v].halfedge = hid;
    }

    #[inline]
    fn set_f_halfedge(&mut self, fid: FaceId, hid: HalfedgeId) {
        self.faces[fid].halfedge = hid;
    }

    #[inline]
    fn set_he_vertex(&mut self, hid: HalfedgeId, vid: VertexId) {
        self.halfedges[hid].vertex = vid;
    }

    #[inline]
    pub fn connect_halfedges(&mut self, hid1: HalfedgeId, hid2: HalfedgeId) {
        self.halfedges[hid1].next = hid2;
        self.halfedges[hid2].prev = hid1;
    }
}

pub trait MeshCore {
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

    fn vertex_mut(&mut self, vid: VertexId) -> &mut Self::Vertex;
    fn halfedge_mut(&mut self, hid: HalfedgeId) -> &mut Self::Halfedge;
    fn face_mut(&mut self, fid: FaceId) -> &mut Self::Face;

    fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2];

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
    fn base(&self) -> &BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A>;
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A>;
}

impl<VP, HP, FP, A: Allocator> HasBaseMesh for BaseMesh<Vertex<VP>, Halfedge<HP>, Face<FP>, A> {
    type A = A;
    type VP = VP;
    type HP = HP;
    type FP = FP;

    #[inline]
    fn base(&self) -> &BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A> {
        self
    }

    #[inline]
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A> {
        self
    }
}

impl<T: HasBaseMesh> MeshCore for T {
    type Vertex = Vertex<T::VP>;
    type Halfedge = Halfedge<T::HP>;
    type Face = Face<T::FP>;

    #[inline]
    fn n_vertices(&self) -> usize {
        self.base().n_vertices
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
    fn n_vertices_capacity(&self) -> usize {
        self.base().vertices.data.len()
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

    fn he_edge(&self, hid: HalfedgeId) -> EdgeId;

    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId;
}
