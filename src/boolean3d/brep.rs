use std::alloc::Allocator;

use crate::{
    geometry::BBox,
    mesh::{Mesh, SurfaceMesh},
    utils::TwoDimArr,
    INVALID_IND,
};

type LoopId = crate::mesh::FaceId;

#[derive(Clone, Copy)]
pub struct FaceId(usize);

pub struct BrepModel<A: Allocator + Copy = std::alloc::Global> {
    mesh: SurfaceMesh<A>,
    loop_faces: Vec<FaceId, A>,
    face_loops: TwoDimArr<LoopId, A>,
    points: Vec<f64, A>,
    pub bbox: BBox,
    face_surfaces: Vec<usize, A>,
}

impl BrepModel {
    pub fn new<U1, T1, U2, T2, T3, T4>(
        loops: T1,
        faces: T2,
        surfaces: T3,
        pts: T4,
        bbox: BBox,
    ) -> Self
    where
        U1: AsRef<[usize]>,
        T1: IntoIterator<Item = U1>,
        U2: AsRef<[usize]>,
        T2: IntoIterator<Item = U2>,
        T3: IntoIterator<Item = usize>,
        T4: IntoIterator<Item = f64>,
    {
        Self::new_in(loops, faces, surfaces, pts, bbox, std::alloc::Global)
    }
}

impl<A: Allocator + Copy> BrepModel<A> {
    pub fn new_in<U1, T1, U2, T2, T3, T4>(
        loops: T1,
        faces: T2,
        surfaces: T3,
        pts: T4,
        bbox: BBox,
        alloc: A,
    ) -> Self
    where
        U1: AsRef<[usize]>,
        T1: IntoIterator<Item = U1>,
        U2: AsRef<[usize]>,
        T2: IntoIterator<Item = U2>,
        T3: IntoIterator<Item = usize>,
        T4: IntoIterator<Item = f64>,
    {
        let mesh = SurfaceMesh::new(loops, alloc);
        let mut face_loops = TwoDimArr::<LoopId, A>::new_in(alloc);
        let mut loop_faces = Vec::new_in(alloc);
        loop_faces.resize(mesh.n_faces(), FaceId(INVALID_IND));
        for (fid, loop_indices) in faces.into_iter().enumerate() {
            for &loop_idx in loop_indices.as_ref() {
                loop_faces[loop_idx] = FaceId(fid);
            }
            face_loops.push(
                loop_indices
                    .as_ref()
                    .iter()
                    .map(|&idx| crate::mesh::FaceId(idx)),
            );
        }
        let mut face_surfaces = Vec::new_in(alloc);
        face_surfaces.extend(surfaces);
        let mut points = Vec::new_in(alloc);
        points.extend(pts);
        Self {
            mesh,
            loop_faces,
            face_loops,
            points,
            bbox,
            face_surfaces,
        }
    }
}
