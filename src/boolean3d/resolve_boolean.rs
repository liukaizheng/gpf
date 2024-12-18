use crate::mesh::{FaceId, SurfaceMesh};

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
    }
}
