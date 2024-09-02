mod adaptive_subdivide;

use std::{alloc::Allocator, collections::HashMap};

use itertools::Itertools;

use crate::{
    geometry::{BBox, Surf, Surface},
    mesh::{FaceId, Mesh, SurfaceMesh, VertexId},
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

struct TetSet {
    mesh: SurfaceMesh,
    tets: Vec<([VertexId; 4], [(FaceId, bool); 4])>,
    points: Vec<f64>,
}

impl TetSet {
    #[inline]
    fn vertices_in<A: Allocator + Copy>(&self, tid: usize, alloc: A) -> Vec<VertexId, A> {
        self.tets[tid].0.to_vec_in(alloc)
    }
}

pub fn boolean3d(first: &SimpleBody, second: &SimpleBody, t: BooleanType) {
    let surfaces = first
        .surfaces
        .iter()
        .chain(second.surfaces.iter())
        .collect_vec();
    let mut bbox = BBox::default();
    bbox.merge(&first.bbox);
    bbox.merge(&second.bbox);
    bbox.scale(1.1);
    let tets = init_mesh(bbox);

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
    let mut tets = Vec::with_capacity(6);
    let mut face_map = HashMap::new();
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
                (fid, true)
            } else {
                let fid: FaceId = triangles.len().into();
                triangles.push(tri);
                face_map.insert(key, fid);
                (fid, false)
            }
        });
        tets.push((t.map(|vid| vid.into()), tet));
    }
    let mesh = SurfaceMesh::new(triangles);
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
    TetSet { mesh, tets, points }
}
