use hashbrown::HashMap;
use tinyvec::TinyVec;

use crate::{boolean3d::{extract_cells::write_chains, write_obj}, mesh::{EdgeId, FaceId, Mesh, SurfaceMesh, VertexId}};

use super::extract_cells::identify_chain_edge;

pub(crate) struct ModelData {
    pub(crate) points: Vec<f64>,
    pub(crate) patches: Vec<Vec<FaceId>>,
    pub(crate) patch_surface_arr: Vec<usize>,
    pub(crate) surface_patches: Vec<Vec<usize>>,
    pub(crate) face_patch_arr: Vec<usize>,
    pub(crate) cells: Vec<Vec<usize>>,
    pub(crate) patch_cell_arr: Vec<usize>,
    pub(crate) mesh: SurfaceMesh,
}

impl ModelData {
    pub(crate) fn resolve(&self) {
        write_obj("123.obj", &self.points, &self.mesh);
        let (non_manifold_vertices, chains, is_chain_edge) = identify_chain_edge(&self.mesh, |fa, fb| {
            let p1 = self.face_patch_arr[fa];
            let p2 = self.face_patch_arr[fb];
            self.patch_surface_arr[p1] == self.patch_surface_arr[p2]
        });

        println!("the number of non_manifold_vertices: {:?}", non_manifold_vertices.len());
        println!("the number of chains: {:?}", chains.len());
        write_chains("chain.obj", &self.points,&self.mesh, &is_chain_edge);
    }
}
