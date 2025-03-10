use std::alloc::Allocator;

use crate::{
    geometry::{BBox, Crv, Surf},
    is_negative,
    mesh::{HoleAwareMesh, ManifoldMesh, Mesh, SurfaceMesh, VertexId},
    point_3, strip_orientation,
    utils::{Bitmask, TwoDimArr},
    INVALID_IND,
};

pub struct UVFace<A: Allocator + Copy> {
    pub(crate) mesh: ManifoldMesh<A>,
    pub(crate) points: Vec<f64, A>,
}

pub struct BrepFace<A: Allocator + Copy> {
    surface_id: usize,
    bbox: BBox,
    uv_face: Option<UVFace<A>>,
}

pub struct NewBrepModel<A: Allocator + Copy = std::alloc::Global> {
    mesh: HoleAwareMesh<A>,
    points: Vec<f64, A>,
    faces: Vec<BrepFace<A>, A>,
    surfaces: Vec<Surf, A>,
    surf_bbox: Vec<BBox, A>,
    edge_curves: Vec<Crv, A>,
    bbox: BBox,
}

impl<A: Allocator + Copy> NewBrepModel<A> {
    pub fn new_in<A1: Allocator + Copy>(
        loops: TwoDimArr<usize, A>,
        face_loops: TwoDimArr<usize, A>,
        points: Vec<f64, A>,
        face_surfaces: &[usize],
        surfaces: Vec<Surf, A>,
        edge_curves: Vec<Crv, A>,
        alloc: A,
    ) {
        let mesh = HoleAwareMesh::new(loops.iter(), face_loops.iter(), alloc);
        for (face, &surface_id) in mesh.faces().zip(face_surfaces) {
            let reversed = is_negative(surface_id);
            let surf = &surfaces[strip_orientation(surface_id)];
        }
    }
}

type LoopId = crate::mesh::FaceId;

#[derive(Clone, Copy)]
pub struct FaceId(usize);

pub struct BrepModel<A: Allocator + Copy = std::alloc::Global> {
    pub(crate) mesh: SurfaceMesh<A>,
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

    pub fn v_mask(&self, vid: VertexId, n_surfaces: usize) -> Bitmask {
        let mut mask = Bitmask::new(n_surfaces);
        for l in self
            .mesh
            .vertex(vid)
            .incoming_halfedges()
            .map(|he| he.face())
        {
            mask.set(strip_orientation(self.face_surfaces[self.loop_faces[*l].0]));
        }
        mask
    }

    #[inline]
    pub fn v_point(&self, vid: VertexId) -> &[f64] {
        point_3(&self.points, vid.0)
    }
}
