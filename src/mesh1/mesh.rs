#![allow(dead_code)]
use std::{alloc::Allocator, marker::PhantomData, ops::{Index, IndexMut}};

use super::element::{ElementId, FaceId, HalfedgeId, VertexId};


pub struct ElementContainer<T, E: ElementId, A: Allocator> {
    data: Vec<T, A>,
    marker: PhantomData<E>,
}

impl<T, E: ElementId, A: Allocator> Index<E> for ElementContainer<T, E, A> {
    type Output = T;

    #[inline]
    fn index(&self, index: E) -> &Self::Output {
        unsafe {
            self.data.get_unchecked(index.index())
        }
    }
}

impl<T, E: ElementId, A: Allocator> IndexMut<E> for ElementContainer<T, E, A> {
    #[inline]
    fn index_mut(&mut self, index: E) -> &mut Self::Output {
        unsafe {
            self.data.get_unchecked_mut(index.index())
        }
    }
}

pub struct BaseMesh<VP, HP, FP, A: Allocator = std::alloc::Global> {
    vertices: ElementContainer<VP, VertexId, A>,
    halfedges: ElementContainer<HP, HalfedgeId, A>,
    faces: ElementContainer<FP, FaceId, A>,
    n_vertices: usize,
    n_halfedges: usize,
    n_faces: usize,
}

impl <VP, HP, FP, A: Allocator + Copy> BaseMesh<VP, HP, FP, A> {
    fn vertex(&self, vid: VertexId) -> &VP {
        &self.vertices[vid]
    }

    fn vertex_mut(&mut self, vid: VertexId) -> &mut VP {
        &mut self.vertices[vid]
    }

    fn halfedge(&self, hid: HalfedgeId) -> &HP {
        &self.halfedges[hid]
    }

    fn halfedge_mut(&mut self, hid: HalfedgeId) -> &mut HP {
        &mut self.halfedges[hid]
    }

    fn face(&self, fid: FaceId) -> &FP {
        &self.faces[fid]
    }

    fn face_mut(&mut self, fid: FaceId) -> &mut FP {
        &mut self.faces[fid]
    }
}

pub struct Vertex<P> {
    halfedge: HalfedgeId,
    property: P,
}

pub struct Halfedge<P> {
    vertex: VertexId,
    next: HalfedgeId,
    prev: VertexId,
    face: FaceId,
    property: P,
}

pub struct Face<P> {
    halfedge: HalfedgeId,
    property: P,
}

impl <VP, HP, FP, A: Allocator + Copy> BaseMesh<Vertex<VP>, Halfedge<HP>, Face<FP>, A> {

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
    fn base_mut(&mut self) -> &mut BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A>;
}

impl <T: HasBaseMesh> MeshCore for T {
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
}
