use hashbrown::{HashMap, hash_map::Entry};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::geometry::Surf;
use crate::math::{square_norm, sub_short};
use crate::point_3;
use crate::{
    boolean3d::{extract_cells::write_chains, write_obj},
    mesh::{FaceId, Mesh, SurfaceMesh, VertexId},
    utils::Bitmask,
};

use super::{BrepModel, extract_cells::identify_chain_edge};

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
    pub(crate) fn resolve(&self, models: Vec<BrepModel>) {
        write_obj("123.obj", &self.points, &self.mesh);
        let (non_manifold_vertices, chains, is_chain_edge) =
            identify_chain_edge(&self.mesh, |fa, fb| {
                let p1 = self.face_patch_arr[fa];
                let p2 = self.face_patch_arr[fb];
                self.patch_surface_arr[p1] == self.patch_surface_arr[p2]
            });

        println!(
            "the number of non_manifold_vertices: {:?}",
            non_manifold_vertices.len()
        );
        println!("the number of chains: {:?}", chains.len());
        write_chains("chain.obj", &self.points, &self.mesh, &is_chain_edge);

        let mut mask_vertices_map = HashMap::<Bitmask, TinyVec<[VertexId; 1]>, _, _>::with_capacity(
            non_manifold_vertices.len(),
        );
        for vid in non_manifold_vertices {
            let mut mask = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
            for f in self
                .mesh
                .vertex(vid)
                .incoming_halfedges()
                .map(|he| he.face())
            {
                mask.set(self.patch_surface_arr[self.face_patch_arr[*f]]);
            }
            match mask_vertices_map.entry(mask) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().push(vid);
                }
                Entry::Vacant(entry) => {
                    entry.insert(TinyVec::from_iter([vid]));
                }
            }
        }

        let model_isomesh_vertices_map = models
            .iter()
            .map(|model| {
                model
                    .mesh
                    .vertices()
                    .map(|m_vert| {
                        let m_vid = *m_vert;
                        let mask = model.v_mask(m_vid, self.surface_patches.len());
                        let mut iso_vid = VertexId::default();
                        let mut min_dist = f64::MAX;
                        let m_pt = model.v_point(m_vid);
                        for &vid in mask_vertices_map.get(&mask).unwrap_or(&TinyVec::new()) {
                            let pt = self.v_point(vid);
                            let dist = square_norm(&sub_short::<3, _>(m_pt, pt));
                            if dist < min_dist {
                                min_dist = dist;
                                iso_vid = vid;
                            }
                        }
                        iso_vid
                    })
                    .collect_vec()
            })
            .collect_vec();

        println!(
            "model_isomesh_vertices_map: {:?}",
            model_isomesh_vertices_map
        );
    }

    fn resolve_face_patches(
        &self,
        surfaces: &[Surf],
        model: &BrepModel,
        non_manifold_vertices: &[VertexId],
    ) {
    }

    #[inline]
    fn v_point(&self, vid: VertexId) -> &[f64] {
        point_3(&self.points, vid.0)
    }
}
