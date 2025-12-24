use crate::{math::{norm, square_norm, sub_short}, mesh1::{BaseEdgeData, HalfedgeData}};

use super::mesh::{BaseVertexData, Mesh};

pub trait VertexPosition<const N: usize> {
    fn position(&self) -> &[f64; N];
}

impl<const N: usize> VertexPosition<N> for [f64; N] {
    #[inline]
    fn position(&self) -> &[f64; N] {
        self
    }
}

impl<const N: usize, P: VertexPosition<N>> VertexPosition<N> for BaseVertexData<P> {
    #[inline]
    fn position(&self) -> &[f64; N] {
        self.property.position()
    }
}

pub trait SquaredEdgeLength {
    fn edge_length_squared(&self) -> f64;
    fn edge_length_squared_mut(&mut self) -> &mut f64;
}

pub trait EdgeLength {
    fn edge_length(&self) -> f64;
    fn edge_length_mut(&mut self) -> &mut f64;
}

impl<P> SquaredEdgeLength for BaseEdgeData<P>
where
    P: SquaredEdgeLength,
{
    #[inline]
    fn edge_length_squared(&self) -> f64 {
        self.property.edge_length_squared()
    }

    #[inline]
    fn edge_length_squared_mut(&mut self) -> &mut f64 {
        self.property.edge_length_squared_mut()
    }
}

impl<P> EdgeLength for super::surface_mesh::BaseEdgeData<P>
where
    P: EdgeLength,
{
    #[inline]
    fn edge_length(&self) -> f64 {
        self.property.edge_length()
    }

    #[inline]
    fn edge_length_mut(&mut self) -> &mut f64 {
        self.property.edge_length_mut()
    }
}

#[inline]
pub fn update_edge_lengths_squared_in_edge_data<const N: usize, M>(mesh: &mut M)
where
    M: Mesh,
    M::VertexData: VertexPosition<N>,
    M::HalfedgeData: HalfedgeData,
    M::EdgeData: SquaredEdgeLength,
{
    for edge in mesh.edges_mut() {
        let [va, vb] = edge.vertices();
        let pa = va.data.position();
        let pb = vb.data.position();
        *edge.data.edge_length_squared_mut() = square_norm(&sub_short::<N, _>(pa, pb));
    }
}

#[inline]
pub fn update_edge_lengths_in_edge_data<const N: usize, M>(mesh: &mut M)
where
    M: Mesh,
    M::VertexData: VertexPosition<N>,
    M::HalfedgeData: HalfedgeData,
    M::EdgeData: EdgeLength,
{
    for edge in mesh.edges_mut() {
        let [va, vb] = edge.vertices();
        let pa = va.data.position();
        let pb = vb.data.position();
        *edge.data.edge_length_mut() = norm(&sub_short::<N, _>(pa, pb));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh1::SurfaceMesh;
    use crate::mesh1::MeshCore;

    #[test]
    fn edge_length_updates_into_edge_property_field() {
        #[derive(Default, Clone)]
        struct EdgeProp {
            square_len: f64,
            len: f64,
        }

        impl SquaredEdgeLength for EdgeProp {
            #[inline]
            fn edge_length_squared(&self) -> f64 {
                self.square_len
            }

            #[inline]
            fn edge_length_squared_mut(&mut self) -> &mut f64 {
                &mut self.square_len
            }
        }

        impl EdgeLength for EdgeProp {
            #[inline]
            fn edge_length(&self) -> f64 {
                self.len
            }

            #[inline]
            fn edge_length_mut(&mut self) -> &mut f64 {
                &mut self.len
            }
        }

        let mut mesh =
            SurfaceMesh::<[f64; 3], (), EdgeProp, (), _>::new_in(vec![vec![0, 1, 2]], std::alloc::Global);

        for v in mesh.vertices_mut() {
            let id = *v.id;
            v.data.property = match id {
                0 => [0.0, 0.0, 0.0],
                1 => [1.0, 0.0, 0.0],
                2 => [0.0, 2.0, 0.0],
                _ => unreachable!(),
            };
        }

        update_edge_lengths_squared_in_edge_data(&mut mesh);
        update_edge_lengths_in_edge_data(&mut mesh);

        for e in mesh.edges() {
            let [va, vb] = mesh.e_vertices(e.id);
            let pa = mesh.vertex_data(va).position();
            let pb = mesh.vertex_data(vb).position();
            let dx = pa[0] - pb[0];
            let dy = pa[1] - pb[1];
            let dz = pa[2] - pb[2];
            let expected = dx * dx + dy * dy + dz * dz;
            let edge_len = expected.sqrt();
            assert!((e.data.property.square_len - expected).abs() < 1e-12);
            assert!((e.data.property.len - edge_len).abs() < 1e-12);
        }
    }
}
