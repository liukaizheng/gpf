mod adaptive_subdivide;
mod ar_in_tet;
mod brep;
mod extract_cells;
mod resolve_boolean;
mod tet_set;

pub use brep::BrepModel;
use brep::SurfRep;
use tinyvec::TinyVec;

use std::{alloc::Allocator, any::TypeId, collections::HashMap};

use adaptive_subdivide::{SurfaceData, adaptive_subdivide};
use ar_in_tet::{Arrangement, IsoVert, extract_iso_surface};
use extract_cells::extract_cells;
use itertools::Itertools;
use tet_set::TetSet;

use crate::geometry::Plane;
use crate::{
    INVALID_IND, Tolerance,
    geometry::{BBox, Surf, UniqueSurface},
    math::dot,
    mesh::{EdgeId, ElementId, FaceId, Mesh, SurfaceMesh, square_edge_length},
    oriented_index, strip_orientation,
};

struct IsoSurfMesh {
    arrangements: Vec<Option<Arrangement>>,
    mesh: SurfaceMesh,
    points: Vec<f64>,
    iso_vertices: Vec<IsoVert>,
    face_positions: Vec<(usize, FaceId)>,
    face_parents: Vec<usize>,
}

pub fn boolean3d<F>(mut models: Vec<BrepModel>, bool_func: F, eps: f64)
where
    F: Fn(&[bool]) -> bool,
{
    let surfaces = merge_same_surfaces(&mut models);
    let surface_bboxes = compute_surface_boxes(&surfaces, &mut models);
    let surface_datum = surface_bboxes
        .into_iter()
        .zip(&surfaces)
        .map(|(sub_bboxes, surf)| SurfaceData {
            surf,
            bbox: BBox::from_boxes(&sub_bboxes),
            sub_bboxes,
        })
        .collect_vec();

    let mut tets = init_mesh(
        BBox::from_boxes(surface_datum.iter().map(|data| &data.bbox)),
        surfaces.len(),
    );
    let vals = adaptive_subdivide(&mut tets, surface_datum, eps * eps);

    let iso_surf_mesh = extract_iso_surface(&tets, vals);
    // write_obj("123.obj", &iso_surf_mesh.points, &iso_surf_mesh.mesh);
    println!("mesh n tets: {}", tets.tet_faces.len());

    let model_data = extract_cells(iso_surf_mesh, &tets, surfaces.len());
    model_data.resolve(models, &surfaces, bool_func);
}

fn init_mesh(bbox: BBox, n_surfaces: usize) -> TetSet {
    const TETS: [[usize; 4]; 6] = [
        [0, 1, 7, 3],
        [7, 0, 5, 1],
        [4, 0, 5, 7],
        [4, 6, 0, 7],
        [0, 7, 6, 2],
        [7, 2, 0, 3],
    ];
    let hash_tri = |mut verts: [usize; 3]| {
        verts.sort();
        (verts[0] << 6) | (verts[1] << 3) | verts[2]
    };
    let mut face_tets = vec![[INVALID_IND; 2]; 18];
    let mut tet_vertices = Vec::with_capacity(6);
    let mut tet_edges = Vec::with_capacity(6);
    let mut tet_faces = Vec::with_capacity(6);
    let mut face_map: HashMap<usize, FaceId> = HashMap::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for t in TETS {
        let tet = [
            [t[1], t[2], t[3]],
            [t[0], t[3], t[2]],
            [t[0], t[1], t[3]],
            [t[0], t[2], t[1]],
        ]
        .map(|tri| {
            let key = hash_tri(tri);
            if let Some(&fid) = face_map.get(&key) {
                face_tets[fid.0][1] = tet_faces.len();
                fid
            } else {
                let fid: FaceId = triangles.len().into();
                face_tets[fid.0][0] = tet_faces.len();
                triangles.push(tri);
                face_map.insert(key, fid);
                fid
            }
        });
        tet_vertices.push(t.map(|vid| vid.into()));
        tet_faces.push(tet);
    }

    let mesh = SurfaceMesh::new(triangles, std::alloc::Global);

    tet_edges.extend(tet_vertices.iter().map(|verts| {
        let mut edges = [EdgeId::default(); 6];
        let mut idx = 0;
        for (&va, &vb) in verts.iter().tuple_combinations() {
            let eid = mesh.e_from_va_vb(va, vb);
            debug_assert!(eid.valid());
            edges[idx] = eid;
            idx += 1;
        }
        debug_assert!(idx == 6);
        edges
    }));

    let points = vec![
        bbox.min[0],
        bbox.min[1],
        bbox.min[2],
        bbox.min[0],
        bbox.min[1],
        bbox.max[2],
        bbox.min[0],
        bbox.max[1],
        bbox.min[2],
        bbox.min[0],
        bbox.max[1],
        bbox.max[2],
        bbox.max[0],
        bbox.min[1],
        bbox.min[2],
        bbox.max[0],
        bbox.min[1],
        bbox.max[2],
        bbox.max[0],
        bbox.max[1],
        bbox.min[2],
        bbox.max[0],
        bbox.max[1],
        bbox.max[2],
    ];
    let square_edge_lengths =
        Vec::from_iter(mesh.edges().map(|e| square_edge_length(&points, *e, &mesh)));
    let surf_indices = vec![TinyVec::from_iter(0..n_surfaces); 4];
    TetSet {
        mesh,
        tet_vertices,
        tet_edges,
        tet_faces,
        face_tets,
        points,
        square_edge_lengths,
        surf_indices,
    }
}

fn merge_same_surfaces<A: Allocator + Copy>(models: &mut [BrepModel<A>]) -> Vec<Surf> {
    let tol = Tolerance::global();
    let mut unique_surface_map = HashMap::<TypeId, UniqueSurface<Surf>>::new();
    let mut surf_positions = Vec::new();
    for (i, model) in models.iter().enumerate() {
        for (j, brep_face) in model.faces.iter().enumerate() {
            let surf = brep_face.surf.as_ori_surf().unwrap().0;
            unique_surface_map
                .entry(surf.real_type_id())
                .or_insert(UniqueSurface::new())
                .add_surf(surf, surf_positions.len(), tol);
            surf_positions.push((i, j));
        }
    }

    let mut surfaces = Vec::new();
    for (surf_type_id, unique_surfaces) in unique_surface_map {
        if surf_type_id == TypeId::of::<Plane>() {
            for (surf, indices) in unique_surfaces.into_surfaces() {
                let sid = surfaces.len();
                let dz = surf.as_plane().unwrap().dz;
                for idx in indices {
                    let (i, j) = surf_positions[idx];
                    let brep_face = &mut models[i].faces[j];
                    let reversed = {
                        let (s, reversed) = brep_face.surf.as_ori_surf().unwrap();
                        (dot(&s.as_plane().unwrap().dz, &dz) < 0.0) ^ reversed
                    };

                    brep_face.surf = SurfRep::I(oriented_index(sid, reversed));
                }
                surfaces.push(surf);
            }
        } else {
            for (surf, indices) in unique_surfaces.into_surfaces() {
                let sid = surfaces.len();
                for idx in indices {
                    let (i, j) = surf_positions[idx];
                    let brep_face = &mut models[i].faces[j];
                    let reversed = brep_face.surf.as_ori_surf().unwrap().1;
                    brep_face.surf = SurfRep::I(oriented_index(sid, reversed));
                }
                surfaces.push(surf);
            }
        }
    }
    surfaces
}

fn compute_surface_boxes<A: Allocator + Copy>(
    surfaces: &[Surf],
    models: &mut [BrepModel<A>],
) -> Vec<TinyVec<[BBox; 1]>> {
    let mut surface_boxes = Vec::with_capacity(surfaces.len());
    surface_boxes.resize(surfaces.len(), TinyVec::<[BBox; 1]>::new());
    for model in models {
        for fid in 0..model.faces.len() {
            let (ori_surf_id, bbox) = model.compute_face_box(fid.into(), surfaces);
            surface_boxes[strip_orientation(ori_surf_id)].push(bbox);
        }
    }
    surface_boxes
}

fn write_obj(name: &str, points: &[f64], mesh: &SurfaceMesh) {
    let mut file = std::fs::File::create(name).unwrap();
    use std::io::Write;
    for i in 0..points.len() / 3 {
        writeln!(
            &mut file,
            "v {} {} {}",
            points[i * 3],
            points[i * 3 + 1],
            points[i * 3 + 2]
        )
        .unwrap();
    }

    for face in mesh.faces() {
        let mut face_str = "f".to_string();
        for v in face.vertices() {
            face_str.push_str(&format!(" {}", v.0 + 1));
        }
        writeln!(&mut file, "{}", face_str).unwrap();
    }
}
