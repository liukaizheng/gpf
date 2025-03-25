use std::alloc::Allocator;

use tinyvec::TinyVec;

use crate::{
    geometry::{BBox, Crv, Surf, Surface},
    is_negative,
    mesh::{ElementIndex, FaceId, HoleAwareMesh, Mesh, VertexId},
    point, strip_orientation,
    utils::{Bitmask, TwoDimArr},
};

pub struct UVFace<A: Allocator + Copy> {
    pub(crate) points: Vec<f64, A>,
    pub(crate) triangles: Vec<usize, A>,
}

impl<A: Allocator + Copy> UVFace<A> {
    pub fn new(points: Vec<f64, A>, triangles: Vec<usize, A>) -> Self {
        Self { points, triangles }
    }
}

pub struct BrepFace<A: Allocator + Copy> {
    surface_id: usize,
    bbox: BBox,
    uv_face: Option<UVFace<A>>,
}

impl<A: Allocator + Copy> BrepFace<A> {
    pub fn new(surface_id: usize, bbox: BBox, uv_face: Option<UVFace<A>>) -> Self {
        Self {
            surface_id,
            bbox,
            uv_face,
        }
    }
}

pub struct BrepModel<A: Allocator + Copy = std::alloc::Global> {
    pub(crate) mesh: HoleAwareMesh<A>,
    points: Vec<f64, A>,
    faces: Vec<BrepFace<A>, A>,
    surfaces: Vec<Surf, A>,
    surf_bboxes: Vec<BBox, A>,
    surf_faces: Vec<TinyVec<[FaceId; 1]>, A>,
    edge_curves: Vec<Crv, A>,
    pub(crate) bbox: BBox,
}

impl<A: Allocator + Copy> BrepModel<A> {
    pub fn new_in<A1: Allocator + Copy>(
        loops: TwoDimArr<usize, A>,
        face_loops: TwoDimArr<usize, A>,
        points: Vec<f64, A>,
        face_surfaces: &[usize],
        surfaces: Vec<Surf, A>,
        edge_curves: Vec<Crv, A>,
        alloc: A,
    ) -> Self {
        let mesh = HoleAwareMesh::new(loops.iter(), face_loops.iter(), alloc);
        let mut brep_faces = Vec::with_capacity_in(mesh.n_faces_capacity(), alloc);
        let mut surf_faces = Vec::with_capacity_in(surfaces.len(), alloc);
        surf_faces.resize(surfaces.len(), TinyVec::new());
        brep_faces.extend(mesh.faces().zip(face_surfaces).map(|(face, &ori_surf_id)| {
            let reversed = is_negative(ori_surf_id);

            let surf_id = strip_orientation(ori_surf_id);
            surf_faces[surf_id].push(*face);
            let surf = &surfaces[surf_id];
            let (face_box, uv_face) = if let Surf::Plane(_p) = surf
                && face.halfedges().all(|he| edge_curves[*he].is_segment())
            {
                (
                    BBox::from_iter(face.vertices().map(|v| point::<3>(&points, v.index()))),
                    None,
                )
            } else {
                let (uv_points, uv_triangles) = if reversed {
                    surf.compute_uv_face(
                        face.wires().map(|wire| {
                            wire.halfedges()
                                .rev()
                                .map(|he| (&edge_curves[*he], he.same_dir()))
                        }),
                        alloc,
                    )
                } else {
                    surf.compute_uv_face(
                        face.wires().map(|wire| {
                            wire.halfedges()
                                .map(|he| (&edge_curves[*he], !he.same_dir()))
                        }),
                        alloc,
                    )
                };
                let uv_face = UVFace::new(uv_points, uv_triangles);
                (
                    surf.compute_box_from_uv_face(&uv_face.points, &uv_face.triangles, alloc),
                    Some(uv_face),
                )
            };
            BrepFace::new(ori_surf_id, face_box, uv_face)
        }));
        let mut surf_bboxes = Vec::with_capacity_in(surfaces.len(), alloc);
        surf_bboxes.extend(surf_faces.iter().map(|faces| {
            BBox::from_boxes(faces.iter().map(|fid: &FaceId| &brep_faces[*fid].bbox))
        }));
        let bbox = BBox::from_boxes(surf_bboxes.iter());
        BrepModel {
            mesh,
            points,
            faces: brep_faces,
            surfaces,
            surf_bboxes,
            surf_faces,
            edge_curves,
            bbox,
        }
    }
    pub fn v_mask(&self, vid: VertexId, n_surfaces: usize) -> Bitmask {
        let mut mask = Bitmask::new(n_surfaces);
        for face in self
            .mesh
            .vertex(vid)
            .incoming_halfedges()
            .map(|he| he.face())
        {
            mask.set(strip_orientation(self.faces[*face].surface_id));
        }
        mask
    }

    #[inline]
    pub fn v_point(&self, vid: VertexId) -> &[f64] {
        point::<3>(&self.points, vid.0)
    }
}
