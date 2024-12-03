use std::{alloc::Allocator, collections::VecDeque};

use bumpalo::Bump;
use hashbrown::{
    hash_map::{DefaultHashBuilder, Entry},
    HashMap, HashSet,
};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    disjoint_set::DisjointSet,
    is_positive,
    mesh::{EdgeId, ElementId, FaceId, Mesh, SurfaceMesh, VertexId},
    oriented_index, point, strip_orientation, twin_index, INVALID_IND,
};

use super::{ar_in_tet::IsoVert, tet_set::TetSet, Arrangement, IsoSurfMesh};

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

fn write_shell(name: &str, iso_surf_mesh: &IsoSurfMesh, shell: &[usize], patches: &[Vec<FaceId>]) {
    let mesh = &iso_surf_mesh.mesh;
    let mut vertex_map = vec![INVALID_IND; mesh.n_vertices_capacity()];

    let mut points = Vec::new();
    for &ori_pid in shell {
        let pid = strip_orientation(ori_pid);
        for &fid in &patches[pid] {
            for v in mesh.face(fid).vertices() {
                let vid = *v;
                if vertex_map[vid] == INVALID_IND {
                    vertex_map[vid] = points.len() / 3;
                    points.extend_from_slice(point(&iso_surf_mesh.points, *vid));
                }
            }
        }
    }

    let mut file = std::fs::File::create(name).unwrap();
    use std::io::Write;

    for p in points.chunks(3) {
        writeln!(&mut file, "v {} {} {}", p[0], p[1], p[2]).unwrap();
    }

    for &ori_pid in shell {
        let pid = strip_orientation(ori_pid);
        for &fid in &patches[pid] {
            let mut v_ids = Vec::new();
            for v in mesh.face(fid).vertices() {
                v_ids.push(vertex_map[*v] + 1);
            }
            if ori_pid & 1 == 0 {
                v_ids.reverse();
            }
            writeln!(&mut file, "f {}", v_ids.iter().join(" ")).unwrap();
        }
    }
}

pub(super) fn extract_components(iso_surf_mesh: IsoSurfMesh, tets: &TetSet) {
    let (chains, is_chain_edge) =
        identify_chain_edge(&iso_surf_mesh.mesh, &iso_surf_mesh.face_parents);
    println!("the n chains is {}", chains.len());
    write_chains("chain.obj", &iso_surf_mesh, &is_chain_edge);

    let (patches, face_patch_arr) = extract_patches(&iso_surf_mesh.mesh, &is_chain_edge);
    println!("the n patches is {}", patches.len());

    let (shells, patch_shell_arr) = extract_shells(
        &iso_surf_mesh,
        tets,
        &chains,
        &patches,
        &face_patch_arr,
        patches.len() << 1,
    );
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
    tets: &TetSet,
    chains: &[Vec<EdgeId>],
    patches: &[Vec<FaceId>],
    face_patch_arr: &[usize],
    n_ori_patches: usize,
) -> (Vec<Vec<usize>>, Vec<usize>) {
    let mut shell_ds = DisjointSet::new(n_ori_patches);
    let mut bump = Bump::new();
    for chain in chains {
        bump.reset();
        order_patches_around_edge(
            iso_surf_mesh,
            tets,
            chain[0],
            face_patch_arr,
            &mut shell_ds,
            &bump,
        );
    }

    // As so far, we don't consider these shells are in different components
    let (shells, patch_shell_arr) = shell_ds.output();
    println!("the n shells is {}", shells.len());

    // disjoin set for components
    let mut comp_ds = DisjointSet::new(shells.len());
    for ori_pa in (0..n_ori_patches).step_by(2) {
        let ori_pb = ori_pa + 1;
        comp_ds.merge(patch_shell_arr[ori_pa], patch_shell_arr[ori_pb]);
    }

    println!("the n comps is {}", comp_ds.n_groups);

    if comp_ds.n_groups > 1 {
        let vert_descent_links = tets.build_descending_vertex_links();
        let (components, shell_to_comp_arr) = comp_ds.output();

        let (tet_vert_iso_vert_arr, component_extremes) =
            find_component_extremes(iso_surf_mesh, tets, &components, patches);
    }

    (shells, patch_shell_arr)
}

struct TetEdgePatchData<A: Allocator + Copy> {
    eid: EdgeId,
    face_to_patches: HashMap<FaceId, TinyVec<[usize; 1]>, DefaultHashBuilder, A>,
}

fn order_patches_around_edge<A: Allocator + Copy>(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    eid: EdgeId,
    face_patches: &[usize],
    ds: &mut DisjointSet,
    alloc: A,
) {
    let mut tet_patch_data = Vec::new_in(alloc);
    let mut tet_to_patch_data_index = HashMap::<usize, usize, _, A>::new_in(alloc);
    for face in iso_surf_mesh.mesh.edge(eid).halfedges().map(|he| he.face()) {
        let fid = *face;
        let (tid, tet_fid) = iso_surf_mesh.face_positions[fid];
        let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();

        let oriented_patch_id = {
            let pid = ar.face_data[tet_fid].pid;
            let oriented_surf_id = *ar.plane_surfaces[pid]
                .iter()
                .find(|&&sid| strip_orientation(sid) == iso_surf_mesh.face_parents[*fid])
                .unwrap();
            oriented_index(face_patches[fid], !is_positive(oriented_surf_id))
        };

        match tet_to_patch_data_index.entry(tid) {
            Entry::Vacant(index) => {
                let [va, vb] = iso_surf_mesh.mesh.e_vertices(eid);
                let tet_eid = ar.mesh.he_edge(ar.find_halfedge(tet_fid, va.0, vb.0));
                debug_assert!(tet_eid.valid());
                let mut vec = TinyVec::new();
                vec.push(oriented_patch_id);
                let mut face_to_patches = HashMap::new_in(alloc);
                face_to_patches.insert(tet_fid, vec);
                index.insert(tet_patch_data.len());
                tet_patch_data.push(TetEdgePatchData {
                    eid: tet_eid,
                    face_to_patches,
                });
            }
            Entry::Occupied(index) => {
                match tet_patch_data[*index.get()].face_to_patches.entry(tet_fid) {
                    Entry::Occupied(mut occupied) => {
                        occupied.get_mut().push(oriented_patch_id);
                    }
                    Entry::Vacant(vacant) => {
                        let mut vec = TinyVec::new();
                        vec.push(oriented_patch_id);
                        vacant.insert(vec);
                    }
                }
            }
        }
    }

    debug_assert!(!tet_to_patch_data_index.is_empty());

    // sort patches for combining unconnected components
    for tet_patches in &mut tet_patch_data {
        for patches in tet_patches.face_to_patches.values_mut() {
            if patches.len() > 1 {
                patches.sort_unstable();
            }
        }
    }

    if tet_to_patch_data_index.len() == 1 {
        for (tid, i) in tet_to_patch_data_index {
            // first and last oriented patch in the tet
            let [first_op, last_op] =
                order_patches_in_tet(iso_surf_mesh, tid, None, &tet_patch_data[i], ds, alloc);
            debug_assert!(first_op != INVALID_IND);
            debug_assert!(first_op != INVALID_IND);
            ds.merge(first_op, last_op);
        }
    } else {
        let tet_eid = {
            let (&tid, &i) = tet_to_patch_data_index.iter().next().unwrap();
            let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
            let [pa, pb] = ar.edges[tet_patch_data[i].eid];
            if pa < 4 && pb < 4 {
                Some(tets.tet_edges[tid][TetSet::tet_edge_index(pa, pb)])
            } else {
                None
            }
        };

        if let Some(tet_eid) = tet_eid {
            let mut first_op = INVALID_IND;
            let mut prev_op = INVALID_IND;
            for [tid, plane_id] in tets.tets_around_edge(tet_eid) {
                if let Some(&idx) = tet_to_patch_data_index.get(&tid) {
                    let tet_edge_and_patches = &tet_patch_data[idx];
                    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();

                    let first_fid = ar
                        .mesh
                        .edge(tet_edge_and_patches.eid)
                        .halfedges()
                        .map(|he| *he.face())
                        .find(|&fid| ar.face_data[fid].pid == plane_id)
                        .unwrap();

                    let [op1, op2] = order_patches_in_tet(
                        &iso_surf_mesh,
                        tid,
                        Some(first_fid),
                        tet_edge_and_patches,
                        ds,
                        alloc,
                    );

                    if prev_op == INVALID_IND {
                        first_op = op1;
                    } else {
                        ds.merge(prev_op, op1);
                    }
                    prev_op = op2;
                }
            }
            debug_assert!(first_op != INVALID_IND);
            debug_assert!(prev_op != INVALID_IND);
            ds.merge(prev_op, first_op);
        } else {
            debug_assert!(tet_to_patch_data_index.len() == 2);
            // edge is on face, there are two related tets
            let edge_tets = {
                let mut iter = tet_to_patch_data_index.iter();
                let tid0 = *iter.next().unwrap().0;
                let tid1 = *iter.next().unwrap().0;
                if *tet_to_patch_data_index.get(&tid0).unwrap() == 0 {
                    [tid0, tid1]
                } else {
                    [tid1, tid0]
                }
            };

            let tet_bdy_faces = [
                (edge_tets[0], tet_patch_data[0].eid),
                (edge_tets[1], tet_patch_data[1].eid),
            ]
            .map(|(tid, eid)| {
                let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                let mesh = &ar.mesh;
                let mut bdy_faces: [FaceId; 2] = [FaceId::default(); 2];
                for (tar, src) in bdy_faces.iter_mut().zip(
                    mesh.edge(eid)
                        .halfedges()
                        .map(|he| *he.face())
                        .filter(|&fid| ar.face_data[fid].pid < 4),
                ) {
                    *tar = src;
                }
                debug_assert!(bdy_faces[0].valid() && bdy_faces[1].valid());
                bdy_faces
            });

            #[derive(PartialEq, Eq)]
            enum Vert {
                /// It's a tet vertex
                Tet(VertexId),
                /// It's a vertex of iso-surface
                ISO(usize),
            }

            let [va, vb] = iso_surf_mesh.mesh.e_vertices(eid);
            let v0 = {
                let tid = edge_tets[0];
                let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                let mesh = &ar.mesh;
                let hid = ar.find_halfedge(tet_bdy_faces[0][0], va.0, vb.0);

                let v = mesh.he_to(mesh.he_next(hid));
                if *v < 4 {
                    Vert::Tet(tets.tet_vertices[tid][*v])
                } else {
                    Vert::ISO(ar.vertices[v].index)
                }
            };
            let v1 = {
                let tid = edge_tets[1];
                let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                let mesh = &ar.mesh;
                let hid = ar.find_halfedge(tet_bdy_faces[1][0], va.0, vb.0);
                let v = mesh.he_from(mesh.he_prev(hid));
                if *v < 4 {
                    Vert::Tet(tets.tet_vertices[tid][*v])
                } else {
                    Vert::ISO(ar.vertices[v].index)
                }
            };

            #[cfg(debug_assertions)]
            {
                let v01 = {
                    let tid = edge_tets[0];
                    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                    let mesh = &ar.mesh;
                    let hid = ar.find_halfedge(tet_bdy_faces[0][1], va.0, vb.0);

                    let v = mesh.he_to(mesh.he_next(hid));
                    if *v < 4 {
                        Vert::Tet(tets.tet_vertices[tid][*v])
                    } else {
                        Vert::ISO(ar.vertices[v].index)
                    }
                };
                let v11 = {
                    let tid = edge_tets[1];
                    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                    let mesh = &ar.mesh;
                    let hid = ar.find_halfedge(tet_bdy_faces[1][1], va.0, vb.0);
                    let v = mesh.he_from(mesh.he_prev(hid));
                    if *v < 4 {
                        Vert::Tet(tets.tet_vertices[tid][*v])
                    } else {
                        Vert::ISO(ar.vertices[v].index)
                    }
                };

                debug_assert!(v0 == v1 || v0 == v11);
                debug_assert!(v01 == v1 || v01 == v11);
            }

            let [pa, pb] = order_patches_in_tet(
                iso_surf_mesh,
                edge_tets[0],
                Some(tet_bdy_faces[0][0]),
                &tet_patch_data[0],
                ds,
                alloc,
            );

            let [pc, pd] = order_patches_in_tet(
                iso_surf_mesh,
                edge_tets[1],
                Some(if v0 == v1 {
                    tet_bdy_faces[1][1]
                } else {
                    tet_bdy_faces[1][0]
                }),
                &tet_patch_data[1],
                ds,
                alloc,
            );

            ds.merge(pb, pc);
            ds.merge(pd, pa);
        }
    }
}

fn order_patches_in_tet<A: Allocator + Copy>(
    iso_surf_mesh: &IsoSurfMesh,
    tid: usize,
    first_fid: Option<FaceId>,
    tet_patches_around_edge: &TetEdgePatchData<A>,
    ds: &mut DisjointSet,
    alloc: A,
) -> [usize; 2] {
    let eid = tet_patches_around_edge.eid;
    let face_to_patches = &tet_patches_around_edge.face_to_patches;
    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
    let mesh = &ar.mesh;
    let mut edge_faces = HashSet::new_in(alloc);
    edge_faces.extend(mesh.edge(eid).halfedges().map(|he| *he.face()));
    let first_fid = if let Some(first_fid) = first_fid {
        first_fid
    } else {
        edge_faces
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
            .unwrap_or(*mesh.edge(eid).halfedge().face())
    };
    let mut curr_fid = first_fid;
    let mut curr_cid = ar.face_data[curr_fid].cells[0];
    debug_assert!(curr_cid != INVALID_IND);

    let mut first_ori_patch = INVALID_IND;

    let mut prev_ori_patch = INVALID_IND;
    loop {
        if curr_fid == first_fid && first_ori_patch != INVALID_IND {
            break;
        }

        if let Some(patches) = face_to_patches.get(&curr_fid) {
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

            if oriented_patches.len() > 1 {
                if oriented_patches[0] != patches[0] {
                    // make sure the first patch (with minimum index as we have sorted the patches) points to the same cell as the face
                    oriented_patches.reverse();
                }
            }

            for (&pa, &pb) in oriented_patches.iter().tuple_windows() {
                ds.merge(twin_index(pa), pb);
            }

            let back_ori_patch = oriented_patches[0];
            let front_ori_patch = twin_index(*oriented_patches.last().unwrap());

            if prev_ori_patch == INVALID_IND {
                first_ori_patch = back_ori_patch;
            } else {
                ds.merge(back_ori_patch, prev_ori_patch);
            }
            prev_ori_patch = front_ori_patch;
        }

        if curr_cid == INVALID_IND {
            break;
        }
        (curr_fid, curr_cid) = {
            let mut next_fid = FaceId::default();
            for &fid in &ar.cell_faces[curr_cid] {
                if fid != curr_fid && edge_faces.contains(&fid) {
                    next_fid = fid;
                    break;
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

    debug_assert!(first_ori_patch != INVALID_IND);
    debug_assert!(prev_ori_patch != INVALID_IND);

    [first_ori_patch, prev_ori_patch]
}

fn find_component_extremes(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    components: &[Vec<usize>],
    patches: &[Vec<FaceId>],
) -> (Vec<VertexId>, Vec<VertexId>) {
    let mut patch_visited = vec![false; patches.len()];
    let mut tet_vert_to_iso_vert_arr = vec![VertexId::default(); tets.mesh.n_vertices()];
    let component_extremes = Vec::from_iter(components.iter().map(|ori_patches| {
        let mut extreme_pt: &[f64] = &[f64::MAX, f64::MAX, f64::MAX];
        let mut extreme_vid = VertexId::default();
        for &op in ori_patches {
            let patch_id = strip_orientation(op);
            if patch_visited[patch_id] {
                continue;
            }
            patch_visited[patch_id] = true;
            for &fid in &patches[patch_id] {
                for v in iso_surf_mesh.mesh.face(fid).vertices() {
                    let vid = *v;
                    match iso_surf_mesh.iso_vertices[vid] {
                        IsoVert::V(tet_vid) => {
                            if tet_vert_to_iso_vert_arr[tet_vid].valid() {
                                continue;
                            }
                            tet_vert_to_iso_vert_arr[tet_vid] = vid;
                            let pt = point(&tets.points, tet_vid.0);
                            if pt.partial_cmp(extreme_pt).unwrap().is_lt() {
                                extreme_pt = pt;
                                extreme_vid = vid;
                            }
                        }
                        IsoVert::ES((tet_eid, _)) => {
                            let [va, vb] = tets.mesh.e_vertices(tet_eid);
                            let p1 = point(&tets.points, va.0);
                            let p2 = point(&tets.points, vb.0);
                            let min_pt = if p1.partial_cmp(&p2).unwrap().is_lt() {
                                p1
                            } else {
                                p2
                            };

                            if min_pt.partial_cmp(extreme_pt).unwrap().is_lt() {
                                extreme_pt = min_pt;
                                extreme_vid = vid;
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
        debug_assert!(extreme_vid.valid());
        extreme_vid
    }));
    (tet_vert_to_iso_vert_arr, component_extremes)
}

fn identify_orient_patches_across_components(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    vert_descent_links: &[VertexId],
    face_patch_arr: &[usize],
    component_extremes: Vec<VertexId>,
) {
    for (comp_id, extreme) in component_extremes.into_iter().enumerate() {
        let(curr_component_patch, next_vid) = match iso_surf_mesh.iso_vertices[extreme] {
            IsoVert::V(tet_vid) => {
                let mut next_vid = vert_descent_links[tet_vid];
                debug_assert!(next_vid.valid());
                (find_component_outer_path_for_vert(
                    iso_surf_mesh,
                    tets,
                    &face_patch_arr,
                    tet_vid,
                    next_vid,
                ), next_vid)
            }
            IsoVert::ES((tet_eid, _)) => {
                let [va, vb] = tets.mesh.e_vertices(tet_eid);
                let p1 = point(&tets.points, va.0);
                let p2 = point(&tets.points, vb.0);
                let min_pt = if p1.partial_cmp(&p2).unwrap().is_lt() {
                    p1
                } else {
                    p2
                };
            }
            _ => panic!("unexpected extreme vertex")
        };
    }
}

fn find_component_outer_path_for_vert(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    face_patch_arr: &[usize],
    start_vid: VertexId,
    end_vid: VertexId,
) -> usize {
    let is_component_tet = |tid: usize| -> bool {
        if let Some(ar) = &iso_surf_mesh.arrangements[tid] {
            let mesh = &ar.mesh;
            for he in mesh.vertex(start_vid).incoming_halfedges() {
                if ar.has_srf_on(*he.face()) {
                    return true;
                }
            }
            false
        } else {
            false
        }
    };

    let start_tid = {
        let eid = tets.mesh.e_from_va_vb(start_vid, end_vid);
        let fid = *tets.mesh.edge(eid).halfedge().face();
        debug_assert!(fid.valid());
        tets.face_tets[fid][0]
    };
    if is_component_tet(start_tid) {
        let t_start_vid_idx = tets.tet_vert_index(start_tid, start_vid);
        let t_end_vid_idx = tets.tet_vert_index(start_tid, end_vid);
        let ar = iso_surf_mesh.arrangements[start_tid].as_ref().unwrap();

        let mut t_bottom_fid = FaceId::default();
        for he in ar.mesh.vertex(t_end_vid_idx.into()).incoming_halfedges() {
            let fid = *he.face();
            if ar.face_data[fid].pid == t_start_vid_idx {
                t_bottom_fid = fid;
                break;
            }
        }
        debug_assert!(t_bottom_fid.valid());

        let cid = ar.face_data[t_bottom_fid].cells[0];
        let t_surf_fid = ar.find_cell_surf_face(cid);
        debug_assert!(t_surf_fid.valid());
        return get_face_patch_in_tet(iso_surf_mesh, ar, t_surf_fid, !ar.is_face_inner_cell(t_surf_fid, cid), face_patch_arr);
    }

    let get_component_tet_and_face = || {
        let mut visited_tets = HashSet::with_capacity(4);
        visited_tets.insert(start_tid);
        let mut queue = VecDeque::new();
        queue.push_back(start_tid);

        while !queue.is_empty() {
            let curr_tid = queue.pop_front().unwrap();
            let vert_pos = tets.tet_vert_index(curr_tid, start_vid);
            let tet_faces = &tets.tet_faces[curr_tid];
            for i in 1..4 {
                let fid = tet_faces[(vert_pos + i) % 4];
                let face_tets = tets.face_tets[fid];
                let tid = face_tets[0] ^ face_tets[1] ^ curr_tid;

                if visited_tets.contains(&tid) {
                    continue;
                }

                if is_component_tet(tid) {
                    return (tid, fid);
                }

                visited_tets.insert(tid);
                queue.push_back(tid);
            }
        }
        (INVALID_IND, FaceId::default())
    };

    let (start_tid, start_fid) = get_component_tet_and_face();
    debug_assert!(start_tid != INVALID_IND);
    let ar = iso_surf_mesh.arrangements[start_tid].as_ref().unwrap();
    let t_start_fid = tets.tet_face_index(start_tid, start_fid).into();
    if ar.has_srf_on(t_start_fid) {
        return get_face_patch_in_tet(iso_surf_mesh, ar, t_start_fid, false, face_patch_arr);
    } else {
        let cid = ar.face_data[t_start_fid].cells[0];
        let t_surf_fid = ar.find_cell_surf_face(cid);
        debug_assert!(t_surf_fid.valid());
        return get_face_patch_in_tet(iso_surf_mesh, ar, t_surf_fid, !ar.is_face_inner_cell(t_surf_fid, cid), face_patch_arr);
    }
}

fn find_component_outer_patch_for_edge(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    face_patch_arr: &[usize],
    t_eid: EdgeId,
    t_curr_vid: VertexId,
    t_next_vid: VertexId,
) -> usize {
    let tid = {
        let fid = *tets.mesh.edge(t_eid).halfedge().face();
        tets.face_tets[fid][0]
    };
    debug_assert!(tid != INVALID_IND);
    debug_assert!(iso_surf_mesh.arrangements[tid].is_some());

    get_patch_from_bottom_face(iso_surf_mesh, tets, tid, face_patch_arr, t_curr_vid, t_next_vid)
}

fn get_patch_from_bottom_face(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    tid: usize,
    face_patch_arr: &[usize],
    t_curr_vid: VertexId,
    t_next_vid: VertexId,
) -> usize {
    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
    let curr_vid_idx = tets.tet_vert_index(tid, t_curr_vid);
    let next_vid = tets.tet_vert_index(tid, t_next_vid).into();

    let mut bottom_fid = FaceId::default();
    for he in ar.mesh.vertex(next_vid).incoming_halfedges() {
        let fid = *he.face();
        if ar.face_data[fid].pid == curr_vid_idx {
            bottom_fid = fid;
            break;
        }
    }
    debug_assert!(bottom_fid.valid());

    let cid = ar.face_data[bottom_fid].cells[0];
    let surf_fid = ar.find_cell_surf_face(cid);
    get_face_patch_in_tet(iso_surf_mesh, ar, surf_fid, !ar.is_face_inner_cell(surf_fid, cid), face_patch_arr)
}

fn get_face_patch_in_tet(
    iso_surf_mesh: &IsoSurfMesh,
    ar: &Arrangement,
    tet_fid: FaceId,
    reversed: bool,
    face_patch_arr: &[usize],
) -> usize {
    let pid = ar.face_data[tet_fid].pid;
    debug_assert!(!ar.plane_surfaces[pid].is_empty());
    let hid = *ar.mesh.face(tet_fid).halfedge();
    let [va, vb] = ar.mesh.he_vertices(hid);
    let vc = ar.mesh.he_to(ar.mesh.he_next(hid));

    let [va, vb, vc] = [va, vb, vc].map(|vid| ar.vertices[vid].index.into());

    let eid = iso_surf_mesh.mesh.e_from_va_vb(va, vb);
    let faces = TinyVec::<[FaceId; 1]>::from_iter(
        iso_surf_mesh.mesh.edge(eid).halfedges().filter_map(|he| {
            if *he.next().to() == vc || *he.prev().from() == vc {
                debug_assert!(iso_surf_mesh.face_positions[*he.face()].1 == tet_fid);
                Some(*he.face())
            } else {
                None
            }
        }),
    );

    debug_assert!(!faces.is_empty());

    let plane_surfs = &ar.plane_surfaces[pid];
    let get_oriented_patch = |fid: FaceId| {
        let surf_id = iso_surf_mesh.face_parents[fid];
        let oriented_sid = *plane_surfs
            .iter()
            .find(|&&sid| strip_orientation(sid) == surf_id)
            .unwrap();
        oriented_index(face_patch_arr[fid], !is_positive(oriented_sid))
    };

    if reversed {
        twin_index(
            faces
                .iter()
                .map(|&fid| get_oriented_patch(fid))
                .max()
                .unwrap(),
        )
    } else {
        faces
            .iter()
            .map(|&fid| get_oriented_patch(fid))
            .min()
            .unwrap()
    }
}
