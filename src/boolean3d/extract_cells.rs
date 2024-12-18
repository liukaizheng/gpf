use std::{alloc::Allocator, collections::VecDeque};

use bumpalo::Bump;
use hashbrown::{
    hash_map::{DefaultHashBuilder, Entry},
    HashMap, HashSet,
};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    boolean3d::tet_set::EDGE_FACE_INDICES,
    is_positive,
    mesh::{EdgeId, ElementId, FaceId, HalfedgeId, Mesh, SurfaceMesh, VertexId},
    oriented_index, point, strip_orientation, twin_index,
    utils::{DisjointSet, TwoDimArr},
    INVALID_IND,
};

use super::{
    ar_in_tet::IsoVert, resolve_boolean::ModelData, tet_set::TetSet, Arrangement, IsoSurfMesh,
};

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

fn write_shell(
    name: &str,
    mesh: &SurfaceMesh,
    all_points: &[f64],
    shell: &[usize],
    patches: &[Vec<FaceId>],
) {
    let mut vertex_map = vec![INVALID_IND; mesh.n_vertices_capacity()];

    let mut points = Vec::new();
    for &ori_pid in shell {
        let pid = strip_orientation(ori_pid);
        for &fid in &patches[pid] {
            for v in mesh.face(fid).vertices() {
                let vid = *v;
                if vertex_map[vid] == INVALID_IND {
                    vertex_map[vid] = points.len() / 3;
                    points.extend_from_slice(point(&all_points, *vid));
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
            for (va, vb) in v_ids[1..].iter().tuple_windows() {
                writeln!(&mut file, "f {} {} {}", v_ids[0], va, vb).unwrap();
            }
        }
    }
}

pub(super) fn extract_cells(
    iso_surf_mesh: IsoSurfMesh,
    tets: &TetSet,
    n_surfaces: usize,
) -> ModelData {
    let (chains, is_chain_edge) =
        identify_chain_edge(&iso_surf_mesh.mesh, &iso_surf_mesh.face_parents);
    println!("the n chains is {}", chains.len());
    // write_chains("chain.obj", &iso_surf_mesh, &is_chain_edge);

    let (patches, face_patch_arr) = extract_patches(&iso_surf_mesh.mesh, &is_chain_edge);
    println!("the n patches is {}", patches.len());

    let (cells, patch_cell_arr) = extract_cells_impl(
        &iso_surf_mesh,
        tets,
        &chains,
        &patches,
        &face_patch_arr,
        patches.len() << 1,
    );
    remove_unused_patches(&iso_surf_mesh, cells, &patches, &patch_cell_arr, n_surfaces)
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
                if edge_visited[curr_eid] {
                    break;
                }
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
fn identify_boundary_patches(mesh: &SurfaceMesh, face_patch_arr: &[usize]) -> Vec<bool> {
    let mut is_boundary_patch_arr = vec![false; face_patch_arr.len()];
    for e in mesh.edges() {
        let he = e.halfedge();
        let next_he = he.sibling();
        if he.eq(&next_he) {
            is_boundary_patch_arr[face_patch_arr[*he.face()]] = true;
        }
    }
    is_boundary_patch_arr
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

fn extract_cells_impl(
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

    let vert_descent_links = tets.build_descending_vertex_links();
    let (components, shell_comp_arr) = comp_ds.output();

    let (tet_vert_iso_elem_arr, component_extremes) =
        find_component_extremes(iso_surf_mesh, tets, &components, &shells, &patches);

    let connect_info = ConnectInfo {
        face_patch_arr: &face_patch_arr,
        patch_shell_arr: &patch_shell_arr,
        shell_component_arr: &shell_comp_arr,
    };

    let outer_patch = get_outer_patch(
        iso_surf_mesh,
        tets,
        &mut shell_ds,
        vert_descent_links,
        connect_info,
        component_extremes,
        tet_vert_iso_elem_arr,
    );

    let (cells, patch_cell_arr) = extract_cells_by_removing_boundary_patches(
        shell_ds,
        &iso_surf_mesh.mesh,
        face_patch_arr,
        n_ori_patches,
        outer_patch,
    );
    println!("the n cells is {}", cells.len());
    // for (shell_id, shell) in cells.iter().enumerate() {
    //     write_shell(
    //         &format!("data/mesh/shell_{}.obj", shell_id),
    //         &iso_surf_mesh.mesh,
    //         &iso_surf_mesh.points,
    //         shell,
    //         patches,
    //     );
    // }

    (cells, patch_cell_arr)
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
                let tet_eid = ar.mesh.he_edge(ar.find_halfedge(tet_fid, va, vb));
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
                ISO(VertexId),
            }

            let [va, vb] = iso_surf_mesh.mesh.e_vertices(eid);
            let v0 = {
                let tid = edge_tets[0];
                let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                let mesh = &ar.mesh;
                let hid = ar.find_halfedge(tet_bdy_faces[0][0], va, vb);

                let v = mesh.he_to(mesh.he_next(hid));
                if *v < 4 {
                    Vert::Tet(tets.tet_vertices[tid][*v])
                } else {
                    Vert::ISO(ar.vertices[v].iso_vid)
                }
            };
            let v1 = {
                let tid = edge_tets[1];
                let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                let mesh = &ar.mesh;
                let hid = ar.find_halfedge(tet_bdy_faces[1][0], va, vb);
                let v = mesh.he_from(mesh.he_prev(hid));
                if *v < 4 {
                    Vert::Tet(tets.tet_vertices[tid][*v])
                } else {
                    Vert::ISO(ar.vertices[v].iso_vid)
                }
            };

            #[cfg(debug_assertions)]
            {
                let v01 = {
                    let tid = edge_tets[0];
                    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                    let mesh = &ar.mesh;
                    let hid = ar.find_halfedge(tet_bdy_faces[0][1], va, vb);

                    let v = mesh.he_to(mesh.he_next(hid));
                    if *v < 4 {
                        Vert::Tet(tets.tet_vertices[tid][*v])
                    } else {
                        Vert::ISO(ar.vertices[v].iso_vid)
                    }
                };
                let v11 = {
                    let tid = edge_tets[1];
                    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
                    let mesh = &ar.mesh;
                    let hid = ar.find_halfedge(tet_bdy_faces[1][1], va, vb);
                    let v = mesh.he_from(mesh.he_prev(hid));
                    if *v < 4 {
                        Vert::Tet(tets.tet_vertices[tid][*v])
                    } else {
                        Vert::ISO(ar.vertices[v].iso_vid)
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

struct ConnectInfo<'a> {
    face_patch_arr: &'a [usize],
    patch_shell_arr: &'a [usize],
    shell_component_arr: &'a [usize],
}

impl<'a> ConnectInfo<'a> {
    #[inline]
    fn get_component_id(&self, fid: FaceId) -> usize {
        self.shell_component_arr[self.patch_shell_arr[self.face_patch_arr[fid] << 1]]
    }
}

#[derive(Clone, Copy)]
enum IsoElem {
    V,
    E(EdgeId),
    None,
}

fn find_component_extremes(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    components: &[Vec<usize>],
    shells: &[Vec<usize>],
    patches: &[Vec<FaceId>],
) -> (Vec<IsoElem>, Vec<VertexId>) {
    let mut patch_visited = vec![false; patches.len()];
    let mut tet_vert_to_iso_elem_arr = vec![IsoElem::None; tets.mesh.n_vertices()];
    let component_extremes = Vec::from_iter(components.iter().map(|comp_shells| {
        let mut extreme_pt: &[f64] = &[f64::MAX, f64::MAX, f64::MAX];
        let mut extreme_vid = VertexId::default();
        for &shell_id in comp_shells {
            for &op in &shells[shell_id] {
                let patch_id = strip_orientation(op);
                if patch_visited[patch_id] {
                    continue;
                }
                patch_visited[patch_id] = true;
                for &fid in &patches[patch_id] {
                    for v in iso_surf_mesh.mesh.face(fid).vertices() {
                        let vid = *v;
                        match iso_surf_mesh.iso_vertices[vid] {
                            IsoVert::V(tet_vid) => match tet_vert_to_iso_elem_arr[tet_vid] {
                                IsoElem::V => {}
                                IsoElem::None | IsoElem::E(_) => {
                                    tet_vert_to_iso_elem_arr[tet_vid] = IsoElem::V;
                                    let pt = point(&tets.points, tet_vid.0);
                                    if pt.partial_cmp(extreme_pt).unwrap().is_lt() {
                                        extreme_pt = pt;
                                        extreme_vid = vid;
                                    }
                                }
                            },
                            IsoVert::ES((tet_eid, _)) => {
                                let [va, vb] = tets.mesh.e_vertices(tet_eid);
                                let p1 = point(&tets.points, va.0);
                                let p2 = point(&tets.points, vb.0);
                                let (min_pt, max_vid) = if p1.partial_cmp(&p2).unwrap().is_lt() {
                                    (p1, vb)
                                } else {
                                    (p2, va)
                                };

                                match tet_vert_to_iso_elem_arr[max_vid] {
                                    IsoElem::V | IsoElem::E(_) => {}
                                    IsoElem::None => {
                                        tet_vert_to_iso_elem_arr[max_vid] = IsoElem::E(tet_eid);
                                    }
                                }

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
        }
        debug_assert!(extreme_vid.valid());
        extreme_vid
    }));
    (tet_vert_to_iso_elem_arr, component_extremes)
}

fn get_outer_patch(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    ds: &mut DisjointSet,
    vert_descent_links: Vec<VertexId>,
    info: ConnectInfo,
    component_extremes: Vec<VertexId>,
    tet_vert_to_iso_elem_arr: Vec<IsoElem>,
) -> usize {
    let mut outer_oriented_patch = INVALID_IND;
    for (comp_id, extreme) in component_extremes.into_iter().enumerate() {
        let (curr_oriented_patch, t_prev_vid, t_next_vid) =
            match iso_surf_mesh.iso_vertices[extreme] {
                IsoVert::V(tet_vid) => {
                    let next_vid = vert_descent_links[tet_vid];
                    debug_assert!(next_vid.valid());
                    (
                        find_component_patch_at_vertex(
                            iso_surf_mesh,
                            tets,
                            Some(ds),
                            &info,
                            comp_id,
                            tet_vid,
                            next_vid,
                        ),
                        tet_vid,
                        next_vid,
                    )
                }
                IsoVert::ES((tet_eid, _)) => {
                    let [va, vb] = tets.mesh.e_vertices(tet_eid);
                    let p1 = point(&tets.points, va.0);
                    let p2 = point(&tets.points, vb.0);
                    let [min_vid, max_vid] = if p1.partial_cmp(&p2).unwrap().is_lt() {
                        [va, vb]
                    } else {
                        [vb, va]
                    };
                    (
                        find_component_patch_at_edge(
                            iso_surf_mesh,
                            tets,
                            Some(ds),
                            &info,
                            comp_id,
                            tet_eid,
                            min_vid,
                        ),
                        max_vid,
                        min_vid,
                    )
                }
                _ => panic!("unexpected extreme vertex"),
            };

        if curr_oriented_patch != INVALID_IND {
            let mut prev_vid = t_prev_vid;
            let mut curr_vid = t_next_vid;
            loop {
                let next_vid = vert_descent_links[curr_vid];
                if !next_vid.valid() {
                    outer_oriented_patch = curr_oriented_patch;
                    break;
                }

                match tet_vert_to_iso_elem_arr[curr_vid] {
                    IsoElem::V => {
                        let next_oriented_patch = find_component_patch_at_vertex(
                            iso_surf_mesh,
                            tets,
                            None,
                            &info,
                            comp_id,
                            curr_vid,
                            prev_vid,
                        );
                        debug_assert!(next_oriented_patch != INVALID_IND);
                        ds.merge(curr_oriented_patch, next_oriented_patch);
                        break;
                    }
                    IsoElem::E(eid) => {
                        let next_oriented_patch = find_component_patch_at_edge(
                            iso_surf_mesh,
                            tets,
                            None,
                            &info,
                            INVALID_IND,
                            eid,
                            curr_vid,
                        );
                        debug_assert!(next_oriented_patch != INVALID_IND);
                        ds.merge(curr_oriented_patch, next_oriented_patch);
                        break;
                    }
                    IsoElem::None => {
                        prev_vid = curr_vid;
                        curr_vid = next_vid;
                    }
                }
            }
        }
    }
    outer_oriented_patch
}

fn find_first_component_intersection_on_edge(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    info: &ConnectInfo,
    t_eid: EdgeId,
    t_start_vid: VertexId,
    comp_id: usize,
) -> (usize, HalfedgeId) {
    let tid = {
        let fid = *tets.mesh.edge(t_eid).halfedge().face();
        debug_assert!(fid.valid());
        tets.face_tets[fid][0]
    };

    debug_assert!(iso_surf_mesh.arrangements[tid].is_some());
    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();
    let mut start_vid = VertexId(tets.tet_vert_index(tid, t_start_vid));

    let edge_index = tets.tet_edges[tid]
        .iter()
        .position(|&eid| t_eid == eid)
        .unwrap();
    let base_edge_planes = EDGE_FACE_INDICES[edge_index];
    let mut prev_eid = EdgeId::default();
    loop {
        let curr_hid = ar
            .mesh
            .vertex(start_vid)
            .outgoing_halfedges()
            .map(|he| *he)
            .find(|&hid| {
                if ar.mesh.he_edge(hid) == prev_eid {
                    return false;
                }
                let eid = ar.mesh.he_edge(hid);
                let mut edge_planes = ar.edges[eid];
                edge_planes.sort();
                base_edge_planes == edge_planes
            });
        if let Some(curr_hid) = curr_hid {
            let vid = ar.mesh.he_to(curr_hid);
            for he in ar.mesh.vertex(vid).incoming_halfedges() {
                let fid = *he.face();
                let iso_fid = ar.face_data[fid].iso_fid;
                if iso_fid.valid() {
                    let curr_comp_id = info.get_component_id(iso_fid);
                    if comp_id == INVALID_IND || curr_comp_id == comp_id {
                        return (tid, curr_hid);
                    }
                }
            }
            start_vid = vid;
            prev_eid = ar.mesh.he_edge(curr_hid);
        } else {
            break;
        }
    }
    debug_assert!(false);
    (tid, HalfedgeId::default())
}

fn find_component_patch_at_vertex(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    ds: Option<&mut DisjointSet>,
    info: &ConnectInfo,
    comp_id: usize,
    t_start_vid: VertexId,
    t_end_vid: VertexId,
) -> usize {
    let is_component_tet = |tid: usize| -> bool {
        if let Some(ar) = &iso_surf_mesh.arrangements[tid] {
            let mesh = &ar.mesh;
            for he in mesh.vertex(t_start_vid).incoming_halfedges() {
                if ar.face_data[*he.face()].iso_fid.valid() {
                    return true;
                }
            }
            false
        } else {
            false
        }
    };

    let start_tid = {
        let eid = tets.mesh.e_from_va_vb(t_start_vid, t_end_vid);
        let fid = *tets.mesh.edge(eid).halfedge().face();
        debug_assert!(fid.valid());
        tets.face_tets[fid][0]
    };
    if is_component_tet(start_tid) {
        let start_vid = VertexId(tets.tet_vert_index(start_tid, t_start_vid));
        let end_vid = VertexId(tets.tet_vert_index(start_tid, t_end_vid));
        debug_assert!(start_vid.valid());
        debug_assert!(end_vid.valid());

        let descent_eid_planes = EDGE_FACE_INDICES[TetSet::tet_edge_index(start_vid.0, end_vid.0)];
        let ar = iso_surf_mesh.arrangements[start_tid].as_ref().unwrap();
        let mut descent_hid = HalfedgeId::default();
        for he in ar.mesh.vertex(start_vid).outgoing_halfedges() {
            let eid = ar.mesh.he_edge(*he);
            let mut edge_planes = ar.edges[eid];
            edge_planes.sort();
            if edge_planes == descent_eid_planes {
                descent_hid = *he;
                break;
            }
        }
        debug_assert!(descent_hid.valid());

        return locate_component_patch_from_halfedge(
            iso_surf_mesh,
            ar,
            ds,
            info,
            comp_id,
            descent_hid,
        );
    }

    let get_component_tet_and_face = || {
        let mut visited_tets = HashSet::with_capacity(4);
        visited_tets.insert(start_tid);
        let mut queue = VecDeque::new();
        queue.push_back(start_tid);

        while !queue.is_empty() {
            let curr_tid = queue.pop_front().unwrap();
            let vert_pos = tets.tet_vert_index(curr_tid, t_start_vid);
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

    let (start_tid, t_start_fid) = get_component_tet_and_face();
    debug_assert!(start_tid != INVALID_IND);
    let ar = iso_surf_mesh.arrangements[start_tid].as_ref().unwrap();
    let start_fid = FaceId(tets.tet_face_index(start_tid, t_start_fid));
    {
        let start_iso_fid = ar.face_data[start_fid].iso_fid;
        if start_iso_fid.valid() {
            let start_comp_id = info.get_component_id(start_iso_fid);
            if (start_comp_id == comp_id) == ds.is_some() {
                return resolve_oriented_face_patch(
                    iso_surf_mesh,
                    ar,
                    start_fid,
                    false,
                    info.face_patch_arr,
                );
            }
        }
    }

    let mut cell_visited = vec![false; ar.cell_faces.len()];
    let mut queue = VecDeque::new();
    queue.push_back(ar.face_data[start_fid].cells[0]);
    cell_visited[*queue.back().unwrap()] = true;
    while !queue.is_empty() {
        let curr_cid = queue.pop_front().unwrap();
        for &fid in &ar.cell_faces[curr_cid] {
            let iso_fid = ar.face_data[fid].iso_fid;
            if iso_fid.valid() {
                let curr_comp_id = info.get_component_id(iso_fid);
                if (curr_comp_id == comp_id) == ds.is_some() {
                    return resolve_oriented_face_patch(
                        iso_surf_mesh,
                        ar,
                        fid,
                        ar.is_face_inner_cell(fid, curr_cid),
                        info.face_patch_arr,
                    );
                }
            }

            for &adj_cid in &ar.face_data[fid].cells {
                if adj_cid != INVALID_IND && !cell_visited[adj_cid] {
                    cell_visited[adj_cid] = true;
                    queue.push_back(adj_cid);
                }
            }
        }
    }
    panic!("can't find the component patch");
}

fn find_component_patch_at_edge(
    iso_surf_mesh: &IsoSurfMesh,
    tets: &TetSet,
    ds: Option<&mut DisjointSet>,
    info: &ConnectInfo,
    comp_id: usize,
    t_eid: EdgeId,
    t_start_vid: VertexId,
) -> usize {
    let (tid, descent_hid) = find_first_component_intersection_on_edge(
        iso_surf_mesh,
        tets,
        &info,
        t_eid,
        t_start_vid,
        comp_id,
    );
    debug_assert!(descent_hid.valid());

    let ar = iso_surf_mesh.arrangements[tid].as_ref().unwrap();

    return locate_component_patch_from_halfedge(iso_surf_mesh, ar, ds, info, comp_id, descent_hid);
}

fn locate_component_patch_from_halfedge(
    iso_surf_mesh: &IsoSurfMesh,
    ar: &Arrangement,
    ds: Option<&mut DisjointSet>,
    info: &ConnectInfo,
    comp_id: usize,
    descent_hid: HalfedgeId,
) -> usize {
    let cid = ar.face_data[ar.mesh.he_face(descent_hid)].cells[0];

    let mut curr_comp_fid = FaceId::default();
    let mut next_comp_fid = FaceId::default();
    for &fid in &ar.cell_faces[cid] {
        if curr_comp_fid.valid() && next_comp_fid.valid() {
            break;
        }
        let iso_fid = ar.face_data[fid].iso_fid;
        if !iso_fid.valid() {
            continue;
        }
        let curr_comp_id = info.get_component_id(iso_fid);
        if curr_comp_id == comp_id {
            curr_comp_fid = fid;
        } else {
            next_comp_fid = fid;
        }
    }
    if let Some(ds) = ds {
        debug_assert!(curr_comp_fid.valid());
        let curr_comp_oriented_patch = resolve_oriented_face_patch(
            iso_surf_mesh,
            ar,
            curr_comp_fid,
            ar.is_face_inner_cell(curr_comp_fid, cid),
            info.face_patch_arr,
        );
        if next_comp_fid.valid() {
            let next_comp_oriented_patch = resolve_oriented_face_patch(
                iso_surf_mesh,
                ar,
                next_comp_fid,
                ar.is_face_inner_cell(next_comp_fid, cid),
                info.face_patch_arr,
            );
            ds.merge(curr_comp_oriented_patch, next_comp_oriented_patch);
            return INVALID_IND;
        } else {
            return curr_comp_oriented_patch;
        }
    } else {
        debug_assert!(next_comp_fid.valid());
        return resolve_oriented_face_patch(
            iso_surf_mesh,
            ar,
            next_comp_fid,
            ar.is_face_inner_cell(next_comp_fid, cid),
            info.face_patch_arr,
        );
    }
}

fn resolve_oriented_face_patch(
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

    let [va, vb, vc] = [va, vb, vc].map(|vid| ar.vertices[vid].iso_vid);

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

fn extract_cells_by_removing_boundary_patches(
    mut ds: DisjointSet,
    mesh: &SurfaceMesh,
    face_patch_arr: &[usize],
    n_ori_patches: usize,
    outer_patch: usize,
) -> (Vec<Vec<usize>>, Vec<usize>) {
    let is_boundary_patch_arr = identify_boundary_patches(mesh, face_patch_arr);
    for pid in 0..(n_ori_patches >> 1) {
        if is_boundary_patch_arr[pid] {
            let oriented_patch = oriented_index(pid, false);
            ds.merge(oriented_patch, outer_patch);
            ds.merge(twin_index(oriented_patch), outer_patch);
        }
    }

    let mut group_to_elements = HashMap::<usize, Vec<usize>>::with_capacity(ds.n_groups - 1);
    let outer_patch_parent = ds.find_set(outer_patch);
    for i in 0..n_ori_patches {
        let p = ds.find_set(i);
        if p != outer_patch_parent {
            group_to_elements.entry(p).or_insert(vec![]).push(i);
        }
    }

    let mut cells = Vec::with_capacity(group_to_elements.len());
    let mut patch_to_cell_arr = vec![INVALID_IND; n_ori_patches];
    for (cid, cell) in group_to_elements.into_values().enumerate() {
        for &element in &cell {
            patch_to_cell_arr[element] = cid;
        }
        cells.push(cell);
    }
    (cells, patch_to_cell_arr)
}

fn remove_unused_patches(
    iso_surf_mesh: &IsoSurfMesh,
    mut cells: Vec<Vec<usize>>,
    old_patches: &[Vec<FaceId>],
    old_patch_cell_arr: &[usize],
    n_surfaces: usize,
) -> ModelData {
    let mut points = Vec::new();
    let mut faces = TwoDimArr::<usize>::new();
    let mut patches = Vec::new();
    let mut face_patch_arr = Vec::new();
    let mut patch_surface_arr = Vec::new();

    let mut point_indices = vec![INVALID_IND; iso_surf_mesh.mesh.n_vertices()];
    let mut face_indices = vec![INVALID_IND; iso_surf_mesh.mesh.n_faces()];
    let mut patch_indices = vec![INVALID_IND; old_patch_cell_arr.len() >> 1];
    for (new_patch_id, old_patch_id) in old_patch_cell_arr
        .chunks(2)
        .enumerate()
        .filter_map(|(patch_id, patch_cells)| {
            if patch_cells[0] != INVALID_IND || patch_cells[1] != INVALID_IND {
                Some(patch_id)
            } else {
                None
            }
        })
        .enumerate()
    {
        patch_indices[old_patch_id] = new_patch_id;
        patch_surface_arr.push(iso_surf_mesh.face_parents[old_patches[old_patch_id][0]]);

        let mut patch = Vec::with_capacity(old_patches[old_patch_id].len());
        for &old_fid in &old_patches[old_patch_id] {
            let new_fid = faces.len();
            face_indices[old_fid] = new_fid;
            patch.push(FaceId(new_fid));
            faces.push(iso_surf_mesh.mesh.face(old_fid).vertices().map(|v| {
                let vid = *v;
                if point_indices[vid] == INVALID_IND {
                    point_indices[vid] = points.len() / 3;
                    points.extend_from_slice(point(&iso_surf_mesh.points, vid.0));
                }
                point_indices[vid]
            }));
            face_patch_arr.push(new_patch_id);
        }
        patches.push(patch);
    }

    let mut surface_patches = vec![Vec::new(); n_surfaces];
    for (patch_id, &surface_id) in patch_surface_arr.iter().enumerate() {
        surface_patches[surface_id].push(patch_id);
    }

    let mesh = SurfaceMesh::new(faces.iter(), std::alloc::Global);
    for cell in cells.iter_mut() {
        for oriented_patch_id in cell.iter_mut() {
            let patch_id = strip_orientation(*oriented_patch_id);
            debug_assert!(patch_indices[patch_id] != INVALID_IND);
            *oriented_patch_id =
                oriented_index(patch_indices[patch_id], !is_positive(*oriented_patch_id));
        }
    }
    let mut patch_cell_arr = vec![INVALID_IND; patches.len() << 1];
    for (cid, cell) in cells.iter().enumerate() {
        for &patch_id in cell {
            patch_cell_arr[patch_id] = cid;
        }
    }
    for (shell_id, shell) in cells.iter().enumerate() {
        write_shell(
            &format!("data/mesh/shell_{}.obj", shell_id),
            &mesh,
            &points,
            shell,
            &patches,
        );
    }

    ModelData {
        points,
        face_patch_arr,
        patch_surface_arr,
        surface_patches,
        patches,
        cells,
        patch_cell_arr,
        mesh,
    }
}
