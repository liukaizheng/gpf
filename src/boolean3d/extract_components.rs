use std::{alloc::Allocator, collections::VecDeque};

use bumpalo::Bump;
use hashbrown::{
    hash_map::{DefaultHashBuilder, Entry},
    HashMap, HashSet,
};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    abs_index,
    disjoint_set::DisjointSet,
    mesh::{EdgeId, ElementId, FaceId, Mesh, SurfaceMesh, VertexId},
    signed_index, twin_index, INVALID_IND,
};

use super::IsoSurfMesh;

fn write_chains(name: &str, iso_surf_mesh: &IsoSurfMesh, is_chain_edge: &[bool]) {
    let file = std::fs::File::create(name).unwrap();
    use std::io::Write;
    for p in iso_surf_mesh.points.chunks(3) {
        writeln!(&file, "v {} {} {}", p[0], p[1], p[2]).unwrap();
    }

    for edge in iso_surf_mesh.mesh.edges() {
        if is_chain_edge[*edge] {
            let he = edge.halfedge();
            let v0 = iso_surf_mesh.mesh.he_from(*he);
            let v1 = iso_surf_mesh.mesh.he_to(*he);
            writeln!(&file, "l {} {}", v0.0 + 1, v1.0 + 1).unwrap();
        }
    }
}

pub(super) fn extract_components(iso_surf_mesh: IsoSurfMesh) {
    let (chains, is_chain_edge) =
        identify_chain_edge(&iso_surf_mesh.mesh, &iso_surf_mesh.face_parents);
    println!("the n chains is {}", chains.len());
    write_chains("chain.obj", &iso_surf_mesh, &is_chain_edge);

    let (patch_faces, face_patches) = extract_patches(&iso_surf_mesh.mesh, &is_chain_edge);
    println!("the n patches is {}", patch_faces.len());

    extract_shells(&iso_surf_mesh, &chains, &face_patches, patch_faces.len());
}

fn identify_chain_edge(
    mesh: &SurfaceMesh,
    face_parents: &[usize],
) -> (Vec<Vec<EdgeId>>, Vec<bool>) {
    let mut is_chain_edge = vec![false; mesh.n_edges_capacity()];
    let mut vertex_edge_map = HashMap::<VertexId, TinyVec<[EdgeId; 2]>>::new();
    for edge in mesh.edges() {
        let he = edge.halfedge();
        let next_he = he.sibling();
        if next_he.ne(&he) {
            if next_he.sibling().ne(&he)
                || face_parents[*he.face()] != face_parents[*next_he.face()]
            {
                let eid = *edge;
                is_chain_edge[eid] = true;
                for vid in mesh.he_vertices(*he) {
                    match vertex_edge_map.entry(vid) {
                        Entry::Occupied(mut entry) => {
                            entry.get_mut().push(eid);
                        }
                        Entry::Vacant(entry) => {
                            let mut vals = TinyVec::new();
                            vals.push(eid);
                            entry.insert(vals);
                        }
                    }
                }
            }
        }
    }

    let propagate_chain =
        |mut curr_vid: VertexId, mut curr_eid: EdgeId, edge_visited: &mut [bool]| {
            let mut chain = vec![curr_eid];
            edge_visited[curr_eid] = true;
            loop {
                let [va, vb] = mesh.e_vertices(curr_eid);
                curr_vid = if va == curr_vid { vb } else { va };

                let candidate_edges = vertex_edge_map.get(&curr_vid).unwrap();
                if candidate_edges.len() != 2 {
                    break;
                }
                curr_eid = if candidate_edges[0] == curr_eid {
                    candidate_edges[1]
                } else {
                    candidate_edges[0]
                };
                chain.push(curr_eid);
                edge_visited[curr_eid] = true;
            }
            chain
        };

    let mut chains = Vec::new();
    let mut edge_visited = vec![false; mesh.n_edges_capacity()];
    for (&vid, edges) in &vertex_edge_map {
        for &eid in edges {
            if edge_visited[eid] {
                continue;
            }
            chains.push(propagate_chain(vid, eid, &mut edge_visited));
        }
    }
    (chains, is_chain_edge)
}

fn extract_patches(mesh: &SurfaceMesh, is_chain_edge: &[bool]) -> (Vec<Vec<FaceId>>, Vec<usize>) {
    let mut patch_faces = Vec::new();
    let mut face_patches = Vec::with_capacity(mesh.n_faces_capacity());
    face_patches.resize(mesh.n_faces_capacity(), INVALID_IND);
    for face in mesh.faces() {
        let fid = *face;
        if face_patches[fid] != INVALID_IND {
            continue;
        }

        let pid = patch_faces.len();
        let mut patch = vec![fid];
        face_patches[fid] = pid;
        let mut queue = VecDeque::new();
        queue.push_back(fid);
        // let _vert = Vec::from_iter(face.vertices().map(|v| *v));

        while !queue.is_empty() {
            let curr_fid = queue.pop_front().unwrap();
            for he in mesh.face(curr_fid).halfedges() {
                let adj_he = he.sibling();
                if adj_he.ne(&he) && !is_chain_edge[*he.edge()] {
                    let adj_fid = *adj_he.face();
                    if face_patches[adj_fid] == INVALID_IND {
                        face_patches[adj_fid] = pid;
                        patch.push(adj_fid);
                        queue.push_back(adj_fid);
                    }
                }
            }
        }
        patch_faces.push(patch);
    }
    (patch_faces, face_patches)
}

fn extract_shells(
    iso_surf_mesh: &IsoSurfMesh,
    chains: &[Vec<EdgeId>],
    face_patches: &[usize],
    n_patches: usize,
) {
    let mut ds = DisjointSet::new(n_patches << 1);
    let mut bump = Bump::new();
    for chain in chains {
        bump.reset();
        order_patches_around_edge(iso_surf_mesh, chain[0], face_patches, &mut ds, &bump);
    }
}

struct TetPatchesAroundEdge<A: Allocator + Copy> {
    eid: EdgeId,
    face_patches_map: HashMap<FaceId, TinyVec<[usize; 1]>, DefaultHashBuilder, A>,
}

fn order_patches_around_edge<A: Allocator + Copy>(
    iso_surf_mesh: &IsoSurfMesh,
    eid: EdgeId,
    face_patches: &[usize],
    ds: &mut DisjointSet,
    alloc: A,
) {
    let mut tet_patches_map = HashMap::<usize, TetPatchesAroundEdge<A>, _, A>::new_in(alloc);
    for face in iso_surf_mesh.mesh.edge(eid).halfedges().map(|he| he.face()) {
        let fid = *face;
        let (tid, tet_fid) = iso_surf_mesh.face_positions[fid];
        let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();

        let signed_patch = {
            let pid = ar.face_data[tet_fid].pid;
            let signed_surf_id = ar.plane_surfaces[pid]
                .iter()
                .find(|&&sid| abs_index(sid) == iso_surf_mesh.face_parents[*fid])
                .unwrap();
            signed_index(face_patches[fid], (signed_surf_id & 1) != 0)
        };

        match tet_patches_map.entry(tid) {
            Entry::Vacant(vacant) => {
                let [va, vb] = iso_surf_mesh.mesh.e_vertices(eid);
                let tet_eid = ar.find_edge(tet_fid, va.0, vb.0);
                debug_assert!(tet_eid.valid());
                let mut vec = TinyVec::new();
                vec.push(signed_patch);
                let mut face_patches_map = HashMap::new_in(alloc);
                face_patches_map.insert(tet_fid, vec);
                vacant.insert(TetPatchesAroundEdge {
                    eid: tet_eid,
                    face_patches_map,
                });
            }
            Entry::Occupied(mut occupied) => {
                match occupied.get_mut().face_patches_map.entry(tet_fid) {
                    Entry::Occupied(mut occupied_entry) => {
                        occupied_entry.get_mut().push(signed_patch);
                    }
                    Entry::Vacant(vacant_entry) => {
                        let mut vec = TinyVec::new();
                        vec.push(signed_patch);
                        vacant_entry.insert(vec);
                    }
                }
            }
        }
    }

    debug_assert!(!tet_patches_map.is_empty());

    if tet_patches_map.len() == 1 {
        for (tid, tet_patches_around_edge) in tet_patches_map {
            order_patches_in_tet(iso_surf_mesh, tid, tet_patches_around_edge, ds, alloc);
        }
    } else {
    }
    println!("there are {} shells", ds.n_groups);
}

fn order_patches_in_tet<A: Allocator + Copy>(
    iso_surf_mesh: &IsoSurfMesh,
    tid: usize,
    tet_patches_around_edge: TetPatchesAroundEdge<A>,
    ds: &mut DisjointSet,
    alloc: A,
) {
    let eid = tet_patches_around_edge.eid;
    let face_patches_map = tet_patches_around_edge.face_patches_map;
    let ar = &iso_surf_mesh.arrangements[tid].as_ref().unwrap();
    let mesh = &ar.mesh;
    let mut edge_faces = HashSet::new_in(alloc);
    edge_faces.extend(mesh.edge(eid).halfedges().map(|he| *he.face()));
    let first_fid = edge_faces
        .iter()
        .find_map(|face| {
            let fid = *face;
            let cells = &ar.face_data[fid].cells;
            if cells[0] == INVALID_IND || cells[1] == INVALID_IND {
                Some(fid)
            } else {
                None
            }
        })
        .unwrap_or(*mesh.edge(eid).halfedge().face());
    let mut curr_fid = first_fid;
    let mut curr_cid = ar.face_data[curr_fid].cells[0];
    debug_assert!(curr_cid != INVALID_IND);

    let mut first_signed_patch = INVALID_IND;

    let mut prev_signed_patch = INVALID_IND;
    loop {
        if curr_cid == INVALID_IND {
            break;
        }

        if curr_fid == first_fid && first_signed_patch != INVALID_IND {
            debug_assert!(prev_signed_patch != INVALID_IND);
            ds.merge(prev_signed_patch, first_signed_patch);
            break;
        }

        if let Some(patches) = face_patches_map.get(&curr_fid) {
            let is_cell_outer_face = {
                let cells = &ar.face_data[curr_fid].cells;
                cells[0] == curr_cid
            };
            let mut oriented_patches = Vec::with_capacity_in(patches.len(), alloc);
            oriented_patches.extend(patches.iter().map(|&patch_id| {
                if is_cell_outer_face {
                    patch_id
                } else {
                    twin_index(patch_id)
                }
            }));

            for (&pa, &pb) in oriented_patches.iter().tuple_windows() {
                ds.merge(twin_index(pa), pb);
            }

            let back_signed_patch = oriented_patches[0];
            let front_signed_patch = twin_index(*oriented_patches.last().unwrap());

            if prev_signed_patch == INVALID_IND {
                first_signed_patch = back_signed_patch;
            } else {
                ds.merge(back_signed_patch, prev_signed_patch);
            }
            prev_signed_patch = front_signed_patch;
        }
        (curr_fid, curr_cid) = {
            let mut next_fid = FaceId::default();
            for &fid in &ar.cell_faces[curr_cid] {
                if fid != curr_fid && edge_faces.contains(&fid) {
                    next_fid = fid;
                }
            }
            debug_assert!(next_fid.valid());
            let cells = &ar.face_data[next_fid].cells;
            debug_assert!(cells[0] == curr_cid || cells[1] == curr_cid);
            (
                next_fid,
                if cells[0] == curr_cid {
                    cells[1]
                } else {
                    cells[0]
                },
            )
        };
    }
}
