mod adaptive_subdivide;
mod tet_set;

use std::{
    alloc::Allocator,
    collections::{HashMap, HashSet},
};

use adaptive_subdivide::adaptive_subdivide;
use itertools::Itertools;
use tet_set::TetSet;

use crate::{
    geometry::{BBox, Surf},
    mesh::{square_edge_length, EdgeId, FaceId, Mesh, SurfaceMesh, VertexId},
    point, INVALID_IND,
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
    let mut tets = init_mesh(bbox);
    adaptive_subdivide(&mut tets, surfaces, eps * eps);

    println!("mesh n  edges: {}", tets.mesh.n_edges());
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
            [t[0], t[1], t[3]],
            [t[1], t[2], t[3]],
            [t[2], t[0], t[3]],
            [t[2], t[1], t[0]],
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

    let mesh = SurfaceMesh::new(triangles);

    tet_edges.extend(tet_faces.iter().map(|faces| {
        let mut set = HashSet::with_capacity(6);
        for &fid in faces {
            for he in mesh.face(fid).halfedges() {
                set.insert(*he.edge());
            }
        }
        debug_assert_eq!(set.len(), 6);
        let mut iter = set.into_iter();
        let mut edges = [EdgeId::default(); 6];
        edges[0] = iter.next().unwrap();
        edges[1] = iter.next().unwrap();
        edges[2] = iter.next().unwrap();
        edges[3] = iter.next().unwrap();
        edges[4] = iter.next().unwrap();
        edges[5] = iter.next().unwrap();
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
    }
}
