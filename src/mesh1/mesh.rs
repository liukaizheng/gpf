use std::{
    alloc::Allocator,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use crate::mesh1::element::{Edge, EdgeHalfedge, EdgeMut};

use super::element::{
    EdgeId, ElementId, Face, FaceData, FaceId, FaceMut, Halfedge, HalfedgeData, HalfedgeId,
    HalfedgeMut, Vertex, VertexData, VertexId, VertexMut,
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

    #[inline]
    pub fn range(&self, start: usize, end: usize) -> impl Iterator<Item = &T> {
        self.data[start..end].iter()
    }

    #[inline]
    pub fn range_mut(&mut self, start: usize, end: usize) -> impl Iterator<Item = &mut T> {
        self.data[start..end].iter_mut()
    }

    #[inline]
    pub fn reserve(&mut self, additional: usize) {
        self.data.reserve(additional);
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
    pub fn vertex_range(&self, start: VertexId, count: usize) -> impl Iterator<Item = &VP> {
        let end = (*start + count).min(self.vertices.len());
        self.vertices.range(*start, end)
    }

    #[inline]
    pub fn vertex_range_mut(
        &mut self,
        start: VertexId,
        count: usize,
    ) -> impl Iterator<Item = &mut VP> {
        let end = (*start + count).min(self.vertices.len());
        self.vertices.range_mut(*start, end)
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
    pub fn halfedge_range(&self, start: HalfedgeId, count: usize) -> impl Iterator<Item = &HP> {
        let end = (*start + count).min(self.halfedges.len());
        self.halfedges.range(*start, end)
    }

    #[inline]
    pub fn halfedge_range_mut(
        &mut self,
        start: HalfedgeId,
        count: usize,
    ) -> impl Iterator<Item = &mut HP> {
        let end = (*start + count).min(self.halfedges.len());
        self.halfedges.range_mut(*start, end)
    }

    #[inline]
    pub fn face(&self, fid: FaceId) -> &FP {
        &self.faces[fid]
    }

    #[inline]
    pub fn face_mut(&mut self, fid: FaceId) -> &mut FP {
        &mut self.faces[fid]
    }

    #[inline]
    pub fn face_range(&self, start: FaceId, count: usize) -> impl Iterator<Item = &FP> {
        let end = (*start + count).min(self.faces.len());
        self.faces.range(*start, end)
    }

    #[inline]
    pub fn face_range_mut(&mut self, start: FaceId, count: usize) -> impl Iterator<Item = &mut FP> {
        let end = (*start + count).min(self.faces.len());
        self.faces.range_mut(*start, end)
    }
}

impl<VP, HP, FP, A: Allocator> BaseMesh<BaseVertexData<VP>, HP, FP, A> {
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

impl<VP: Default + Clone, HP: Default, FP, A: Allocator> BaseMesh<VP, HP, FP, A> {
    #[inline]
    pub fn new_vertices(&mut self, n: usize) -> VertexId {
        let ret = VertexId(self.vertices.len().into());
        self.vertices.data.resize(*ret + n, VP::default());
        self.n_vertices += n;
        ret
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
pub struct BaseVertexData<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> VertexData for BaseVertexData<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }

    #[inline]
    fn set_halfedge(&mut self, hid: HalfedgeId) {
        self.halfedge = hid;
    }
}

impl<P> BaseVertexData<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        BaseVertexData { halfedge, property }
    }

    #[inline]
    fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

#[derive(Default, Clone)]
pub struct BaseHalfedgeData<P> {
    pub vertex: VertexId,
    pub next: HalfedgeId,
    pub prev: HalfedgeId,
    pub face: FaceId,
    pub property: P,
}

impl<P> BaseHalfedgeData<P> {
    #[inline]
    pub fn new(
        vertex: VertexId,
        next: HalfedgeId,
        prev: HalfedgeId,
        face: FaceId,
        property: P,
    ) -> Self {
        BaseHalfedgeData {
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

impl<P> HalfedgeData for BaseHalfedgeData<P> {
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
pub struct BaseFaceData<P> {
    pub halfedge: HalfedgeId,
    pub property: P,
}

impl<P> BaseFaceData<P> {
    #[inline]
    pub fn new(halfedge: HalfedgeId, property: P) -> Self {
        BaseFaceData { halfedge, property }
    }

    #[inline]
    pub fn valid(&self) -> bool {
        self.halfedge.valid()
    }
}

impl<P> FaceData for BaseFaceData<P> {
    #[inline]
    fn halfedge(&self) -> HalfedgeId {
        self.halfedge
    }

    #[inline]
    fn set_halfedge(&mut self, halfedge: HalfedgeId) {
        self.halfedge = halfedge;
    }
}

impl<V: VertexData, H: HalfedgeData, F: FaceData, A: Allocator> BaseMesh<V, H, F, A> {
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
    pub fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.halfedge(hid).face()
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
    type VertexData;
    type HalfedgeData;
    type FaceData;

    fn n_vertices(&self) -> usize;
    fn n_halfedges(&self) -> usize;
    fn n_faces(&self) -> usize;

    fn n_vertices_capacity(&self) -> usize;
    fn n_halfedges_capacity(&self) -> usize;
    fn n_faces_capacity(&self) -> usize;

    fn vertex_reserve(&mut self, additional: usize);
    fn halfedge_reserve(&mut self, additional: usize);
    fn face_reserve(&mut self, additional: usize);

    fn vertex_data(&self, vid: VertexId) -> &Self::VertexData;
    fn halfedge_data(&self, hid: HalfedgeId) -> &Self::HalfedgeData;
    fn face_data(&self, fid: FaceId) -> &Self::FaceData;

    fn vertex_datum(&self) -> impl Iterator<Item = &Self::VertexData>;
    fn halfedge_datum(&self) -> impl Iterator<Item = &Self::HalfedgeData>;
    fn face_datum(&self) -> impl Iterator<Item = &Self::FaceData>;

    fn vertex_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::VertexData>;
    fn halfedge_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::HalfedgeData>;
    fn face_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::FaceData>;

    fn vertex_data_mut(&mut self, vid: VertexId) -> &mut Self::VertexData;
    fn halfedge_data_mut(&mut self, hid: HalfedgeId) -> &mut Self::HalfedgeData;
    fn face_data_mut(&mut self, fid: FaceId) -> &mut Self::FaceData;

    fn vertex_range(
        &'_ self,
        vid: VertexId,
        count: usize,
    ) -> impl Iterator<Item = Vertex<'_, Self>>;
    fn halfedge_range(
        &'_ self,
        hid: HalfedgeId,
        count: usize,
    ) -> impl Iterator<Item = Halfedge<'_, Self>>;
    fn face_range(&'_ self, fid: FaceId, count: usize) -> impl Iterator<Item = Face<'_, Self>>;

    fn vertex_range_mut(
        &'_ mut self,
        vid: VertexId,
        count: usize,
    ) -> impl Iterator<Item = VertexMut<'_, Self>>;
    fn halfedge_range_mut(
        &'_ mut self,
        hid: HalfedgeId,
        count: usize,
    ) -> impl Iterator<Item = HalfedgeMut<'_, Self>>;
    fn face_range_mut(
        &'_ mut self,
        fid: FaceId,
        count: usize,
    ) -> impl Iterator<Item = FaceMut<'_, Self>>;

    fn vertex(&'_ self, vid: VertexId) -> Vertex<'_, Self>;
    fn halfedge(&'_ self, hid: HalfedgeId) -> Halfedge<'_, Self>;
    fn face(&'_ self, fid: FaceId) -> Face<'_, Self>;

    fn vertex_mut(&'_ mut self, vid: VertexId) -> VertexMut<'_, Self>;
    fn halfedge_mut(&'_ mut self, hid: HalfedgeId) -> HalfedgeMut<'_, Self>;
    fn face_mut(&'_ mut self, fid: FaceId) -> FaceMut<'_, Self>;

    fn vertices(&'_ self) -> impl Iterator<Item = Vertex<'_, Self>>;
    fn halfedges(&'_ self) -> impl Iterator<Item = Halfedge<'_, Self>>;
    fn faces(&'_ self) -> impl Iterator<Item = Face<'_, Self>>;

    fn vertices_mut(&'_ mut self) -> impl Iterator<Item = VertexMut<'_, Self>>;
    fn halfedges_mut(&'_ mut self) -> impl Iterator<Item = HalfedgeMut<'_, Self>>;
    fn faces_mut(&'_ mut self) -> impl Iterator<Item = FaceMut<'_, Self>>;

    fn he_from(&self, hid: HalfedgeId) -> VertexId;
    fn he_to(&self, hid: HalfedgeId) -> VertexId;
    fn he_vertices(&self, hid: HalfedgeId) -> [VertexId; 2];
    fn he_prev(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_next(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_face(&self, hid: HalfedgeId) -> FaceId;

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
    ) -> &BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    >;
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    >;
}

impl<VP, HP, FP, A: Allocator> HasBaseMesh
    for BaseMesh<BaseVertexData<VP>, BaseHalfedgeData<HP>, BaseFaceData<FP>, A>
{
    type A = A;
    type VP = VP;
    type HP = HP;
    type FP = FP;

    #[inline]
    fn base(
        &self,
    ) -> &BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    > {
        self
    }

    #[inline]
    fn base_mut(
        &mut self,
    ) -> &mut BaseMesh<
        BaseVertexData<Self::VP>,
        BaseHalfedgeData<Self::HP>,
        BaseFaceData<Self::FP>,
        Self::A,
    > {
        self
    }
}

impl<T: HasBaseMesh> MeshCore for T {
    type VertexData = BaseVertexData<T::VP>;
    type HalfedgeData = BaseHalfedgeData<T::HP>;
    type FaceData = BaseFaceData<T::FP>;

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
    fn vertex_reserve(&mut self, additional: usize) {
        self.base_mut().vertices.reserve(additional);
    }

    #[inline]
    fn halfedge_reserve(&mut self, additional: usize) {
        self.base_mut().halfedges.reserve(additional);
    }

    #[inline]
    fn face_reserve(&mut self, additional: usize) {
        self.base_mut().faces.reserve(additional);
    }

    #[inline]
    fn vertex_data(&self, vid: VertexId) -> &Self::VertexData {
        self.base().vertex(vid)
    }

    #[inline]
    fn halfedge_data(&self, hid: HalfedgeId) -> &Self::HalfedgeData {
        self.base().halfedge(hid)
    }

    #[inline]
    fn face_data(&self, fid: FaceId) -> &Self::FaceData {
        self.base().face(fid)
    }

    #[inline]
    fn vertex_data_mut(&mut self, vid: VertexId) -> &mut Self::VertexData {
        self.base_mut().vertex_mut(vid)
    }

    #[inline]
    fn halfedge_data_mut(&mut self, hid: HalfedgeId) -> &mut Self::HalfedgeData {
        self.base_mut().halfedge_mut(hid)
    }

    #[inline]
    fn face_data_mut(&mut self, fid: FaceId) -> &mut Self::FaceData {
        self.base_mut().face_mut(fid)
    }

    #[inline]
    fn vertex_datum(&self) -> impl Iterator<Item = &Self::VertexData> {
        self.base().vertices.iter()
    }

    #[inline]
    fn halfedge_datum(&self) -> impl Iterator<Item = &Self::HalfedgeData> {
        self.base().halfedges.iter()
    }

    #[inline]
    fn face_datum(&self) -> impl Iterator<Item = &Self::FaceData> {
        self.base().faces.iter()
    }

    #[inline]
    fn vertex_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::VertexData> {
        self.base_mut().vertices.iter_mut()
    }

    #[inline]
    fn halfedge_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::HalfedgeData> {
        self.base_mut().halfedges.iter_mut()
    }

    #[inline]
    fn face_datum_mut(&mut self) -> impl Iterator<Item = &mut Self::FaceData> {
        self.base_mut().faces.iter_mut()
    }

    #[inline]
    fn vertex_range(
        &'_ self,
        vid: VertexId,
        count: usize,
    ) -> impl Iterator<Item = Vertex<'_, Self>> {
        self.base()
            .vertex_range(vid, count)
            .zip(*vid..*vid + count)
            .map(|(data, vid)| Vertex::new_with_data(vid.into(), data, self))
    }

    #[inline]
    fn halfedge_range(
        &'_ self,
        hid: HalfedgeId,
        count: usize,
    ) -> impl Iterator<Item = Halfedge<'_, Self>> {
        self.base()
            .halfedge_range(hid, count)
            .zip(*hid..*hid + count)
            .map(|(data, vid)| Halfedge::new_with_data(vid.into(), data, self))
    }

    #[inline]
    fn face_range(&'_ self, fid: FaceId, count: usize) -> impl Iterator<Item = Face<'_, Self>> {
        self.base()
            .face_range(fid, count)
            .zip(*fid..*fid + count)
            .map(|(data, vid)| Face::new_with_data(vid.into(), data, self))
    }

    #[inline]
    fn vertex_range_mut(
        &'_ mut self,
        vid: VertexId,
        count: usize,
    ) -> impl Iterator<Item = VertexMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        self.base_mut()
            .vertex_range_mut(vid, count)
            .zip(*vid..*vid + count)
            .map(move |(data, vid)| unsafe {
                VertexMut::new_with_data(vid.into(), data, &mut *mesh_ptr)
            })
    }

    #[inline]
    fn halfedge_range_mut(
        &'_ mut self,
        hid: HalfedgeId,
        count: usize,
    ) -> impl Iterator<Item = HalfedgeMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        self.base_mut()
            .halfedge_range_mut(hid, count)
            .zip(*hid..*hid + count)
            .map(move |(data, hid)| unsafe {
                HalfedgeMut::new_with_data(hid.into(), data, &mut *mesh_ptr)
            })
    }

    #[inline]
    fn face_range_mut(
        &'_ mut self,
        fid: FaceId,
        count: usize,
    ) -> impl Iterator<Item = FaceMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        self.base_mut()
            .face_range_mut(fid, count)
            .zip(*fid..*fid + count)
            .map(move |(data, fid)| unsafe {
                FaceMut::new_with_data(fid.into(), data, &mut *mesh_ptr)
            })
    }

    #[inline]
    fn vertex(&'_ self, vid: VertexId) -> Vertex<'_, Self> {
        Vertex::new(vid, self)
    }

    #[inline]
    fn halfedge(&'_ self, hid: HalfedgeId) -> Halfedge<'_, Self> {
        Halfedge::new(hid, self)
    }

    #[inline]
    fn face(&'_ self, fid: FaceId) -> Face<'_, Self> {
        Face::new(fid, self)
    }

    #[inline]
    fn vertex_mut(&'_ mut self, vid: VertexId) -> VertexMut<'_, Self> {
        VertexMut::new(vid, self)
    }

    #[inline]
    fn halfedge_mut(&'_ mut self, hid: HalfedgeId) -> HalfedgeMut<'_, Self> {
        HalfedgeMut::new(hid, self)
    }

    #[inline]
    fn face_mut(&'_ mut self, fid: FaceId) -> FaceMut<'_, Self> {
        FaceMut::new(fid, self)
    }

    #[inline]
    fn vertices(&'_ self) -> impl Iterator<Item = Vertex<'_, Self>> {
        self.vertex_datum()
            .zip(0..self.n_vertices_capacity())
            .filter_map(|(data, vid)| {
                if data.halfedge.valid() {
                    Some(Vertex::new_with_data(vid.into(), data, self))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn halfedges(&'_ self) -> impl Iterator<Item = Halfedge<'_, Self>> {
        self.halfedge_datum()
            .zip(0..self.n_halfedges_capacity())
            .filter_map(|(data, hid)| {
                if data.vertex.valid() {
                    Some(Halfedge::new_with_data(hid.into(), data, self))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn faces(&'_ self) -> impl Iterator<Item = Face<'_, Self>> {
        self.face_datum()
            .zip(0..self.n_faces_capacity())
            .filter_map(|(data, fid)| {
                if data.halfedge.valid() {
                    Some(Face::new_with_data(fid.into(), data, self))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn vertices_mut(&'_ mut self) -> impl Iterator<Item = VertexMut<'_, Self>> {
        let n_vertices_capacity = self.n_vertices_capacity();
        let mesh_ptr = self as *mut Self;
        self.vertex_datum_mut()
            .zip(0..n_vertices_capacity)
            .filter_map(move |(data, vid)| unsafe {
                if data.halfedge.valid() {
                    Some(VertexMut::new_with_data(vid.into(), data, &mut *mesh_ptr))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn halfedges_mut(&'_ mut self) -> impl Iterator<Item = HalfedgeMut<'_, Self>> {
        let n_halfedges_capacity = self.n_halfedges_capacity();
        let mesh_ptr = self as *mut Self;
        self.halfedge_datum_mut()
            .zip(0..n_halfedges_capacity)
            .filter_map(move |(data, hid)| unsafe {
                if data.vertex.valid() {
                    Some(HalfedgeMut::new_with_data(hid.into(), data, &mut *mesh_ptr))
                } else {
                    None
                }
            })
    }

    #[inline]
    fn faces_mut(&'_ mut self) -> impl Iterator<Item = FaceMut<'_, Self>> {
        let mesh_ptr = self as *mut Self;
        let n_faces_capacity = self.n_faces_capacity();
        self.face_datum_mut()
            .zip(0..n_faces_capacity)
            .filter_map(move |(data, fid)| unsafe {
                if data.halfedge.valid() {
                    Some(FaceMut::new_with_data(fid.into(), data, &mut *mesh_ptr))
                } else {
                    None
                }
            })
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
    fn he_face(&self, hid: HalfedgeId) -> FaceId {
        self.base().he_face(hid)
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
    type EdgeData;

    fn n_edges(&self) -> usize;
    fn n_edges_capacity(&self) -> usize;
    fn edge_reserve(&mut self, additional: usize);

    fn edge_datum(&'_ self) -> impl Iterator<Item = &Self::EdgeData>;
    fn edge_datum_mut(&'_ mut self) -> impl Iterator<Item = &mut Self::EdgeData>;

    fn edge_data(&self, eid: EdgeId) -> &Self::EdgeData;
    fn edge_data_mut(&mut self, eid: EdgeId) -> &mut Self::EdgeData;

    fn edge_range(&'_ self, eid: EdgeId, count: usize) -> impl Iterator<Item = Edge<'_, Self>>;

    fn edge_range_mut(
        &'_ mut self,
        eid: EdgeId,
        count: usize,
    ) -> impl Iterator<Item = EdgeMut<'_, Self>>;

    fn edge(&'_ self, eid: EdgeId) -> Edge<'_, Self>;
    fn edge_mut(&'_ mut self, eid: EdgeId) -> EdgeMut<'_, Self>;

    fn edges(&'_ self) -> impl Iterator<Item = Edge<'_, Self>>;
    fn edges_mut(&'_ mut self) -> impl Iterator<Item = EdgeMut<'_, Self>>;

    fn he_edge(&self, hid: HalfedgeId) -> EdgeId;

    fn he_sibling(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_incoming_next(&self, hid: HalfedgeId) -> HalfedgeId;
    fn he_from_oppo_vertex(&self, fid: FaceId, vid: VertexId) -> HalfedgeId;

    fn e_halfedge(&self, eid: EdgeId) -> HalfedgeId;
    fn e_from_vertices(&self, va: VertexId, vb: VertexId) -> EdgeId;

    fn f_halfedge(&self, fid: FaceId) -> HalfedgeId;
}

pub(super) fn edge_from_vertices<M: Mesh>(mesh: &M, va: VertexId, vb: VertexId) -> EdgeId
where
    M::VertexData: VertexData,
    M::HalfedgeData: HalfedgeData,
{
    for edge in mesh.vertex(va).edges() {
        let he = edge.halfedge();
        let v1 = he.to();
        if v1.id == va {
            if he.from().id == vb {
                return edge.id;
            }
        } else if v1.id == vb {
            if he.from().id == va {
                return edge.id;
            }
        }
    }
    EdgeId::default()
}

pub(super) fn halfedge_from_oppo_vertex<M: Mesh>(mesh: &M, fid: FaceId, vid: VertexId) -> HalfedgeId
where
    M::VertexData: VertexData,
    M::HalfedgeData: HalfedgeData,
    M::FaceData: FaceData,
{
    let face = mesh.face(fid);
    let he = face.halfedges().find(|he| he.to().id == vid).unwrap();
    he.prev().id
}
