use super::{element::ElementId, mesh::{BaseVertexData, Mesh, MeshCore}};

pub trait VertexPosition3 {
    fn position(&self) -> &[f64; 3];
}

impl VertexPosition3 for [f64; 3] {
    #[inline]
    fn position(&self) -> &[f64; 3] {
        self
    }
}

impl<P: VertexPosition3> VertexPosition3 for BaseVertexData<P> {
    #[inline]
    fn position(&self) -> &[f64; 3] {
        self.property.position()
    }
}

pub trait EdgeLengthSquaredField {
    fn edge_length_squared(&self) -> f64;
    fn set_edge_length_squared(&mut self, value: f64);
}

pub trait EdgeLengthField {
    fn edge_length(&self) -> f64;
    fn set_edge_length(&mut self, value: f64);
}

impl EdgeLengthSquaredField for f64 {
    #[inline]
    fn edge_length_squared(&self) -> f64 {
        *self
    }

    #[inline]
    fn set_edge_length_squared(&mut self, value: f64) {
        *self = value;
    }
}

impl EdgeLengthField for f64 {
    #[inline]
    fn edge_length(&self) -> f64 {
        *self
    }

    #[inline]
    fn set_edge_length(&mut self, value: f64) {
        *self = value;
    }
}

impl<P> EdgeLengthSquaredField for super::surface_mesh::BaseEdgeData<P>
where
    P: EdgeLengthSquaredField,
{
    #[inline]
    fn edge_length_squared(&self) -> f64 {
        self.property.edge_length_squared()
    }

    #[inline]
    fn set_edge_length_squared(&mut self, value: f64) {
        self.property.set_edge_length_squared(value);
    }
}

impl<P> EdgeLengthField for super::surface_mesh::BaseEdgeData<P>
where
    P: EdgeLengthField,
{
    #[inline]
    fn edge_length(&self) -> f64 {
        self.property.edge_length()
    }

    #[inline]
    fn set_edge_length(&mut self, value: f64) {
        self.property.set_edge_length(value);
    }
}

#[inline]
pub fn update_edge_lengths_squared_in_edge_data<M>(mesh: &mut M)
where
    M: Mesh,
    M::VertexData: VertexPosition3,
    M::EdgeData: EdgeLengthSquaredField,
{
    let mut values = vec![0.0; mesh.n_edges_capacity()];
    for edge in mesh.edges() {
        let eid = edge.id;
        let [va, vb] = mesh.e_vertices(eid);
        if !va.valid() || !vb.valid() {
            continue;
        }
        let pa = mesh.vertex_data(va).position();
        let pb = mesh.vertex_data(vb).position();
        let dx = pa[0] - pb[0];
        let dy = pa[1] - pb[1];
        let dz = pa[2] - pb[2];
        values[eid.index()] = dx * dx + dy * dy + dz * dz;
    }

    for edge in mesh.edges_mut() {
        let eid = edge.id;
        edge.data.set_edge_length_squared(values[eid.index()]);
    }
}

#[inline]
pub fn update_edge_lengths_in_edge_data<M>(mesh: &mut M)
where
    M: Mesh,
    M::VertexData: VertexPosition3,
    M::EdgeData: EdgeLengthField,
{
    let mut values = vec![0.0; mesh.n_edges_capacity()];
    for edge in mesh.edges() {
        let eid = edge.id;
        let [va, vb] = mesh.e_vertices(eid);
        if !va.valid() || !vb.valid() {
            continue;
        }
        let pa = mesh.vertex_data(va).position();
        let pb = mesh.vertex_data(vb).position();
        let dx = pa[0] - pb[0];
        let dy = pa[1] - pb[1];
        let dz = pa[2] - pb[2];
        values[eid.index()] = (dx * dx + dy * dy + dz * dz).sqrt();
    }

    for edge in mesh.edges_mut() {
        let eid = edge.id;
        edge.data.set_edge_length(values[eid.index()]);
    }
}

pub trait EdgeLengthInDataExt: Mesh
where
    Self::VertexData: VertexPosition3,
{
    #[inline]
    fn update_edge_lengths_squared_in_edge_data(&mut self)
    where
        Self::EdgeData: EdgeLengthSquaredField,
    {
        update_edge_lengths_squared_in_edge_data(self)
    }

    #[inline]
    fn update_edge_lengths_in_edge_data(&mut self)
    where
        Self::EdgeData: EdgeLengthField,
    {
        update_edge_lengths_in_edge_data(self)
    }
}

impl<M> EdgeLengthInDataExt for M
where
    M: Mesh,
    M::VertexData: VertexPosition3,
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh1::SurfaceMesh;

    #[test]
    fn edge_length_updates_into_edge_property_field() {
        #[derive(Default, Clone)]
        struct EdgeProp {
            square_len: f64,
            len: f64,
        }

        impl EdgeLengthSquaredField for EdgeProp {
            #[inline]
            fn edge_length_squared(&self) -> f64 {
                self.square_len
            }

            #[inline]
            fn set_edge_length_squared(&mut self, value: f64) {
                self.square_len = value;
            }
        }

        impl EdgeLengthField for EdgeProp {
            #[inline]
            fn edge_length(&self) -> f64 {
                self.len
            }

            #[inline]
            fn set_edge_length(&mut self, value: f64) {
                self.len = value;
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
