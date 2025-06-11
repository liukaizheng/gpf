use std::alloc::Allocator;

use crate::{
    geometry::{BBox, Crv, Curve, Surf, Surface, segment::Segment},
    is_negative,
    mesh::{EdgeId, ElementIndex, FaceId, HalfedgeId, HoleAwareMesh, Mesh, VertexId},
    point, strip_orientation,
    utils::Bitmask,
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

pub enum SurfRep {
    S((Surf, bool)),
    I(usize),
}

impl SurfRep {
    #[inline]
    pub fn as_index(&self) -> Option<usize> {
        match self {
            SurfRep::I(idx) => Some(*idx),
            _ => None,
        }
    }

    #[inline]
    pub fn as_ori_surf(&self) -> Option<(&Surf, bool)> {
        match self {
            SurfRep::S((surf, reversed)) => Some((surf, *reversed)),
            _ => None,
        }
    }
}

pub struct BrepFace<A: Allocator + Copy> {
    pub(crate) surf: SurfRep,
    uv_face: Option<UVFace<A>>,
}

impl<A: Allocator + Copy> BrepFace<A> {
    pub fn new(surface: Surf, reversed: bool) -> Self {
        Self {
            surf: SurfRep::S((surface, reversed)),
            uv_face: None,
        }
    }
}

fn get_halfedge_curve<'a, A: Allocator + Copy>(
    mesh: &HoleAwareMesh<A>,
    edge_curves: &'a [Crv],
    hid: HalfedgeId,
    surf_same_dir: bool,
) -> (&'a Crv, bool) {
    let crv = &edge_curves[mesh.he_edge(hid)];
    (crv, mesh.he_same_dir(hid) ^ surf_same_dir)
}

pub struct BrepModel<A: Allocator + Copy = std::alloc::Global> {
    pub(crate) mesh: HoleAwareMesh<A>,
    points: Vec<f64, A>,
    pub(crate) faces: Vec<BrepFace<A>, A>,
    pub(crate) edge_curves: Vec<Crv, A>,
    alloc: A,
}

impl<A: Allocator + Copy> BrepModel<A> {
    pub fn new_in<
        U1: AsRef<[usize]>,
        T1: IntoIterator<Item = U1>,
        U2: AsRef<[usize]>,
        T2: IntoIterator<Item = U2>,
    >(
        loops: T1,
        face_loops: T2,
        points: Vec<f64, A>,
        face_surfaces: Vec<(Surf, bool), A>,
        two_verts_curves: Vec<([usize; 2], Crv), A>,
        alloc: A,
    ) -> Self {
        let mesh = HoleAwareMesh::new(loops.into_iter(), face_loops.into_iter(), alloc);
        let mut brep_faces = Vec::with_capacity_in(mesh.n_faces_capacity(), alloc);
        let mut edge_curves = Vec::with_capacity_in(mesh.n_edges_capacity(), alloc);
        edge_curves.resize(mesh.n_edges_capacity(), Crv::default());
        for ([va, vb], curve) in two_verts_curves {
            let [va, vb] = [va.into(), vb.into()];
            let eid = mesh.e_from_va_vb(va, vb);
            if *mesh.edge(eid).halfedge().from() == va {
                edge_curves[eid] = curve;
            } else {
                edge_curves[eid] = curve.reversed();
            }
        }
        for edge in mesh.edges() {
            let eid = *edge;
            if edge_curves[eid].is_none() {
                let [va, vb] = mesh.e_vertices(eid);
                let s = point::<3>(&points, va.index());
                let e = point::<3>(&points, vb.index());
                edge_curves[eid] =
                    Crv::Segment(Segment::new([s[0], s[1], s[2]], [e[0], e[1], e[2]]));
            }
        }

        brep_faces.extend(
            face_surfaces
                .into_iter()
                .map(|(surf, reversed)| BrepFace::new(surf, reversed)),
        );

        BrepModel {
            mesh,
            points,
            faces: brep_faces,
            edge_curves,
            alloc,
        }
    }

    pub fn compute_face_box(&mut self, fid: FaceId, surfaces: &[Surf]) -> (usize, BBox) {
        let brep_face = &mut self.faces[fid];
        let ori_surf_id = brep_face.surf.as_index().unwrap();
        let reversed = is_negative(ori_surf_id);

        let surf_id = strip_orientation(ori_surf_id);
        let surf = &surfaces[surf_id];
        let face = self.mesh.face(fid.into());
        let (face_box, uv_face) = if let Surf::Plane(plane) = surf {
            (
                plane.compute_box(
                    face.wires()
                        .next()
                        .unwrap()
                        .halfedges()
                        .map(|he| &self.edge_curves[*he.edge()]),
                    self.alloc,
                ),
                None,
            )
        } else {
            let (uv_points, uv_triangles) = if reversed {
                surf.compute_uv_face(
                    face.wires().map(|wire| {
                        wire.halfedges()
                            .rev()
                            .map(|he| get_halfedge_curve(&self.mesh, &self.edge_curves, *he, false))
                    }),
                    self.alloc,
                )
            } else {
                surf.compute_uv_face(
                    face.wires().map(|wire| {
                        wire.halfedges()
                            .map(|he| get_halfedge_curve(&self.mesh, &self.edge_curves, *he, true))
                    }),
                    self.alloc,
                )
            };
            let uv_face = UVFace::new(uv_points, uv_triangles);
            (
                surf.compute_box_from_uv_face(
                    &uv_face.points,
                    &uv_face.triangles,
                    std::alloc::Global,
                ),
                Some(uv_face),
            )
        };
        brep_face.uv_face = uv_face;
        (ori_surf_id, face_box)
    }

    pub fn vert_mask(&self, vid: VertexId, n_surfaces: usize) -> Bitmask {
        let mut mask = Bitmask::new(n_surfaces);
        for face in self
            .mesh
            .vertex(vid)
            .incoming_halfedges()
            .map(|he| he.face())
        {
            mask.set(strip_orientation(
                self.faces[*face].surf.as_index().unwrap(),
            ));
        }
        mask
    }

    pub fn edge_mask(&self, eid: EdgeId, n_surfaces: usize) -> Bitmask {
        let mut mask = Bitmask::new(n_surfaces);
        for face in self.mesh.edge(eid).halfedges().map(|he| he.face()) {
            mask.set(strip_orientation(
                self.faces[*face].surf.as_index().unwrap(),
            ));
        }
        mask
    }

    #[inline]
    pub fn oriented_face_surface(&self, fid: FaceId) -> usize {
        self.faces[fid].surf.as_index().unwrap()
    }

    #[inline]
    pub fn v_point(&self, vid: VertexId) -> &[f64] {
        point::<3>(&self.points, vid.0)
    }
}
