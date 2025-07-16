use std::alloc::Allocator;

use super::element::{EdgeId, HalfedgeId};

use super::mesh::{BaseMesh, ElementContainer, Face, Halfedge, HasBaseMesh, Mesh, Vertex};

struct HEdge<P> {
    edge: EdgeId,
    property: P,
}

struct Edge<P> {
    halfedge: HalfedgeId,
    property: P,
}

pub struct SurfaceMesh<VP, HP, EP, FP, A: Allocator> {
    base_mesh: BaseMesh<Vertex<VP>, Halfedge<HEdge<HP>>, Face<FP>, A>,
    edges: ElementContainer<Edge<EP>, EdgeId, A>,
    n_edges: usize,
}

impl<VP, HP, EP, FP, A: Allocator> HasBaseMesh
    for SurfaceMesh<VP, HP, EP, FP, A>
{
    type A = A;
    type VP = VP;
    type HP = HEdge<HP>;
    type FP = FP;

    fn base(&self) -> &BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A> {
        &self.base_mesh
    }

    fn base_mut(&mut self) -> &mut BaseMesh<Vertex<Self::VP>, Halfedge<Self::HP>, Face<Self::FP>, Self::A> {
        todo!()
    }

}
