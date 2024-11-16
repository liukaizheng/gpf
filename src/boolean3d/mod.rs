mod adaptive_subdivide;
mod ar_in_tet;
mod tet_set;

use std::collections::HashMap;

use adaptive_subdivide::adaptive_subdivide;
use ar_in_tet::{extract_iso_surface, Arrangement, InterPt};
use itertools::Itertools;
use tet_set::TetSet;

use crate::{
    geometry::{BBox, Surf},
    mesh::{square_edge_length, EdgeId, ElementId, FaceId, Mesh, SurfaceMesh},
    INVALID_IND,
};

pub struct SimpleBody {
    surfaces: Vec<Surf>,
    bbox: BBox,
}

impl SimpleBody {
    pub fn new(surfaces: Vec<Surf>, bbox: BBox) -> SimpleBody {
        SimpleBody { surfaces, bbox }
    }
}

pub enum BooleanType {
    Union,
    Intersection,
    Difference,
}

struct IsoSurfMesh {
    arrangements: Vec<Option<Arrangement>>,
    mesh: SurfaceMesh,
    points: Vec<f64>,
    iso_vertices: Vec<InterPt>,
    face_positions: Vec<(usize, FaceId)>,
    face_parents: Vec<usize>,
}

pub fn boolean3d(first: &SimpleBody, second: &SimpleBody, t: BooleanType, eps: f64) {
    let surfaces = first
        .surfaces
        .iter()
        .chain(second.surfaces.iter())
        .collect_vec();
    let mut bbox = BBox::default();
    bbox.merge(&first.bbox);
    bbox.merge(&second.bbox);
    bbox.scale(1.1);
    bbox.min = [-0.5, -0.5, -0.5];
    bbox.max = [2.0, 2.0, 2.0];
    let mut tets = init_mesh(bbox);
    let vals = adaptive_subdivide(&mut tets, surfaces, eps * eps);

    let iso_surf_mesh = extract_iso_surface(&tets, vals);
    write_obj("123.obj", &iso_surf_mesh.points, &iso_surf_mesh.mesh);

    println!("mesh n tets: {}", tets.tet_faces.len());
}

fn init_mesh(bbox: BBox) -> TetSet {
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
    TetSet {
        mesh,
        tet_vertices,
        tet_edges,
        tet_faces,
        face_tets,
        points,
        square_edge_lengths,
        tet_indices: vec![0; 6],
    }
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
