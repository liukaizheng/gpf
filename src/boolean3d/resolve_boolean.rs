use std::cell::RefCell;

use hashbrown::{HashMap, HashSet};
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::boolean3d::extract_cells::write_shell;
use crate::geometry::{Surf, Surface};
use crate::graphcut::{ArcBuilder, MaxFlow, PushRelabelFifo};
use crate::math::{cross, norm, square_norm, sub_short};
use crate::mesh::{EdgeId, ElementId, HalfedgeId, HoleAwareMesh};
use crate::utils::TwoDimArr;
use crate::{
    INVALID_IND, face_area_2d, is_negative, is_positive, oriented_index, point, point_3,
    strip_orientation, twin_index,
};
use crate::{
    boolean3d::{extract_cells::write_chains, write_obj},
    mesh::{FaceId, Mesh, SurfaceMesh, VertexId},
    utils::Bitmask,
};

use super::{BrepModel, extract_cells::identify_chain_edge};

pub(crate) struct ModelData {
    points: Vec<f64>,
    patches: Vec<Vec<FaceId>>,
    patch_surface_arr: Vec<usize>,
    surface_patches: Vec<Vec<usize>>,
    face_patch_arr: Vec<usize>,
    cells: Vec<Vec<usize>>,
    patch_cell_arr: Vec<usize>,
    mesh: SurfaceMesh,
    patch_areas: RefCell<Option<Vec<f64>>>,
}

impl ModelData {
    pub(crate) fn new(
        points: Vec<f64>,
        patches: Vec<Vec<FaceId>>,
        patch_surface_arr: Vec<usize>,
        surface_patches: Vec<Vec<usize>>,
        face_patch_arr: Vec<usize>,
        cells: Vec<Vec<usize>>,
        patch_cell_arr: Vec<usize>,
        mesh: SurfaceMesh,
    ) -> Self {
        Self {
            points,
            patches,
            patch_surface_arr,
            surface_patches,
            face_patch_arr,
            cells,
            patch_cell_arr,
            mesh,
            patch_areas: RefCell::new(None),
        }
    }

    pub(crate) fn resolve<F>(&self, models: Vec<BrepModel>, surfaces: &[Surf], func: F)
    where
        F: Fn(&[bool]) -> bool,
    {
        write_obj("123.obj", &self.points, &self.mesh);
        let (non_manifold_vertices, chains, edge_chain_indices) =
            identify_chain_edge(&self.mesh, |fa, fb| {
                self.face_patch_arr[fa] == self.face_patch_arr[fb]
            });
        let (patch_mesh, chain_masks, patch_edge_chain_indices) =
            self.build_patch_mesh(&chains, surfaces, &non_manifold_vertices);

        println!(
            "the number of non_manifold_vertices: {:?}",
            non_manifold_vertices.len()
        );
        println!("the number of chains: {:?}", chains.len());
        write_chains("chain.obj", &self.points, &self.mesh, &edge_chain_indices);

        let mut mask_vertices_map =
            HashMap::<Bitmask, TinyVec<[VertexId; 1]>>::with_capacity(patch_mesh.n_vertices());

        let mut add_into_mask_vertices_map = |mask: Bitmask, vid: VertexId| {
            mask_vertices_map.entry(mask).or_default().push(vid);
        };

        let mut chain_data =
            ChainData::new(chains, patch_edge_chain_indices, &self.mesh, &self.points);

        for v in patch_mesh.vertices() {
            let mut mask = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
            for he in v.incoming_halfedges() {
                mask.set(self.patch_surface_arr[*he.face()]);
            }
            if mask.n_elements() <= 3 {
                add_into_mask_vertices_map(mask, *v);
            } else {
                for k in 0..mask.n_elements() {
                    let mut m = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
                    for elements in mask.iter_set_bits().combinations(k) {
                        m.set_from_iter(elements);
                    }
                    add_into_mask_vertices_map(m, *v);
                }
            }
        }

        /* Initialize boundary patche masks,
        2 means unvisited
        0 means boundary oriented patch has same direction with patch
        1 means boundary oriented halfedge has opposite direction with patch */
        let mut bdy_patch_masks = vec![2u8; self.patches.len()];
        let mut cell_visited = vec![false; self.cells.len()];
        let model_cells = models
            .iter()
            .map(|model| {
                self.identify_model_cells(
                    model,
                    &patch_mesh,
                    &chain_masks,
                    &mask_vertices_map,
                    &non_manifold_vertices,
                    &mut bdy_patch_masks,
                    &mut cell_visited,
                    &mut chain_data,
                )
            })
            .collect_vec();

        println!("model_cells: {:?}", model_cells);

        let result_oriented_patches = self.get_result_oriented_patches(model_cells, func);
        write_shell(
            "result.obj",
            &self.mesh,
            &self.points,
            &result_oriented_patches,
            &self.patches,
        );
    }

    fn identify_model_cells(
        &self,
        model: &BrepModel,
        patch_mesh: &HoleAwareMesh<std::alloc::Global>,
        chain_masks: &[Bitmask],
        mask_vertices_map: &HashMap<Bitmask, TinyVec<[VertexId; 1]>>,
        non_manifold_vertices: &[VertexId],
        bdy_patch_masks: &mut [u8],
        cell_visited: &mut [bool],
        chain_data: &mut ChainData,
    ) -> Vec<usize> {
        let iso_vertices = model
            .mesh
            .vertices()
            .map(|m_vert| {
                let m_vid = *m_vert;
                let mask = model.vert_mask(m_vid, self.surface_patches.len());
                let mut iso_vid = VertexId::default();
                let mut min_dist = f64::MAX;
                let m_pt = model.v_point(m_vid);
                for &vid in mask_vertices_map.get(&mask).unwrap_or(&TinyVec::new()) {
                    let pt = self.v_point(non_manifold_vertices[vid]);
                    let dist = square_norm(&sub_short::<3, _>(m_pt, pt));
                    if dist < min_dist {
                        min_dist = dist;
                        iso_vid = vid;
                    }
                }
                iso_vid
            })
            .collect_vec();
        let edge_chains_arr = Vec::from_iter(model.mesh.edges().map(|edge| {
            let [va, vb] = model.mesh.e_vertices(*edge).map(|v| iso_vertices[v]);
            if !va.valid() || !vb.valid() {
                None
            } else {
                let edge_mask = model.edge_mask(*edge, self.surface_patches.len());
                let mut edge_chains = Vec::new();
                if self.propagate_edge_chain(
                    &edge_mask,
                    va,
                    vb,
                    patch_mesh,
                    chain_masks,
                    &mut edge_chains,
                ) {
                    Some(edge_chains)
                } else {
                    None
                }
            }
        }));

        let get_ori_patch_indices = |model_fid, occupied_patches: Vec<FaceId>| {
            Vec::from_iter(occupied_patches.into_iter().map(|fid| {
                oriented_index(*fid, is_negative(model.oriented_face_surface(model_fid)))
            }))
        };

        /* Initialize boundary edge masks,
        2 means unvisited
        0 means boundary halfedge has same direction with edge
        1 means boundary halfedge has opposite direction with edge */
        let mut bdy_edge_masks = vec![2u8; patch_mesh.n_edges()];
        let mut patch_face_visited = vec![false; patch_mesh.n_faces()];
        let model_face_ori_patches = model
            .mesh
            .faces()
            .map(|model_face| {
                let face_ori_surf_id = model.oriented_face_surface(*model_face);
                let face_surf_id = strip_orientation(face_ori_surf_id);
                if model_face
                    .halfedges()
                    .all(|he| edge_chains_arr[*he.edge()].is_some())
                {
                    let mut start_patch_fid = FaceId::default();
                    for model_he in model_face.halfedges() {
                        let model_he_same_dir = model_he.same_dir();
                        for &hid in edge_chains_arr[*model_he.edge()].as_ref().unwrap() {
                            let eid = patch_mesh.he_edge(hid);
                            if patch_mesh.he_same_dir(hid) == model_he_same_dir {
                                bdy_edge_masks[eid] = 0;
                            } else {
                                bdy_edge_masks[eid] = 1;
                            }
                            if !start_patch_fid.valid() {
                                for he in patch_mesh.edge(eid).halfedges() {
                                    if patch_mesh.he_same_dir(*he) == (bdy_edge_masks[eid] == 0) {
                                        let patch_fid = *he.face();
                                        if self.patch_surface_arr[patch_fid] == face_surf_id {
                                            start_patch_fid = patch_fid;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let occupied_patches = self.find_face_occupied_patches(
                        patch_mesh,
                        start_patch_fid,
                        &mut bdy_edge_masks,
                        &mut patch_face_visited,
                    );

                    // reset mask
                    for model_he in model_face.halfedges() {
                        for &hid in edge_chains_arr[*model_he.edge()].as_ref().unwrap() {
                            bdy_edge_masks[patch_mesh.he_edge(hid)] = 2;
                        }
                    }
                    if let Some(occupied_patches) = occupied_patches {
                        return get_ori_patch_indices(*model_face, occupied_patches);
                    }
                };
                get_ori_patch_indices(
                    *model_face,
                    self.find_face_occupied_patches_fallback(
                        face_surf_id,
                        patch_mesh,
                        &model.mesh,
                        &edge_chains_arr,
                        model_face.halfedges().map(|he| *he),
                        chain_data,
                    ),
                )
            })
            .collect_vec();

        println!("model face oriented patches {:?}", model_face_ori_patches);
        for &ori_pid in model_face_ori_patches.iter().flatten() {
            let pid = strip_orientation(ori_pid);
            if is_positive(ori_pid) {
                bdy_patch_masks[pid] = 1;
            } else {
                bdy_patch_masks[pid] = 0;
            }
        }

        // the patches of model are oriented opposite to its interior,
        // the patches of cell are oriented to its interior
        let model_cells = self.find_model_occupied_cells(
            self.patch_cell_arr[twin_index(model_face_ori_patches[0][0])],
            bdy_patch_masks,
            cell_visited,
        );
        for &ori_pid in model_face_ori_patches.iter().flatten() {
            bdy_patch_masks[strip_orientation(ori_pid)] = 2;
        }

        if let Some(model_cells) = model_cells {
            model_cells
        } else {
            self.find_model_occupied_cells_fallback(model_face_ori_patches.into_iter().flatten())
        }
    }

    fn find_face_occupied_patches(
        &self,
        patch_mesh: &HoleAwareMesh<std::alloc::Global>,
        start_fid: FaceId,
        bdy_edge_masks: &mut [u8],
        patch_face_visited: &mut [bool],
    ) -> Option<Vec<FaceId>> {
        let mut occupied_patches = vec![start_fid];
        patch_face_visited[start_fid] = true;
        let mut idx = 0;
        let mut valid = true;
        loop {
            if idx >= occupied_patches.len() {
                break;
            }
            let curr_fid = occupied_patches[idx];
            idx += 1;
            for he in patch_mesh.face(curr_fid).halfedges() {
                let mask = bdy_edge_masks[*he.edge()];
                if mask == 2 {
                    let twin_fid = he
                        .edge()
                        .halfedges()
                        .filter(|h| **h != *he)
                        .map(|h| *h.face())
                        .find(|&f| self.patch_surface_arr[f] == self.patch_surface_arr[curr_fid]);
                    if let Some(twin_fid) = twin_fid
                        && !patch_face_visited[twin_fid]
                    {
                        patch_face_visited[twin_fid] = true;
                        occupied_patches.push(twin_fid);
                    }
                } else if (mask == 1) == he.same_dir() {
                    valid = false;
                    break;
                }
            }
            if !valid {
                break;
            }
        }
        for i in 0..occupied_patches.len() {
            patch_face_visited[occupied_patches[i]] = false;
        }
        if !valid { None } else { Some(occupied_patches) }
    }

    /// use graphcut to label patches
    fn find_face_occupied_patches_fallback<T: IntoIterator<Item = HalfedgeId>, M: Mesh>(
        &self,
        surf_id: usize,
        patch_mesh: &M,
        model_mesh: &M,
        model_edge_chains_arr: &[Option<Vec<HalfedgeId>>],
        model_face_halfedges: T,
        chain_data: &mut ChainData,
    ) -> Vec<FaceId> {
        let patch_to_idx_map = HashMap::<FaceId, usize>::from_iter(
            self.surface_patches[surf_id]
                .iter()
                .enumerate()
                .map(|(idx, &pid)| (pid.into(), idx)),
        );

        let mut edge_is_boundary_map =
            HashMap::<EdgeId, bool>::with_capacity(patch_to_idx_map.len() * 2);

        let mut edge_len_sum = 0.0;
        for &fid in &self.surface_patches[surf_id] {
            for he in patch_mesh.face(fid.into()).halfedges() {
                let eid = *he.edge();
                if !edge_is_boundary_map.contains_key(&eid) {
                    edge_is_boundary_map.insert(eid, false);
                    edge_len_sum += chain_data.get_edge_length(eid);
                }
            }
        }

        let mut internal_costs = vec![0.0; patch_to_idx_map.len()];
        let mut external_costs = vec![0.0; patch_to_idx_map.len()];
        for model_hid in model_face_halfedges {
            let model_eid = model_mesh.he_edge(model_hid);
            let model_he_same_dir = model_mesh.he_same_dir(model_hid);
            if let Some(chains) = &model_edge_chains_arr[model_eid] {
                for &chain_hid in chains {
                    let chain_eid = patch_mesh.he_edge(chain_hid);
                    edge_is_boundary_map.insert(chain_eid, true);
                    for he in patch_mesh.edge(chain_eid).halfedges() {
                        let hid = *he;
                        let patch_fid = *he.face();
                        if let Some(&idx) = patch_to_idx_map.get(&patch_fid) {
                            if patch_mesh.hes_same_dir(chain_hid, hid) == model_he_same_dir {
                                external_costs[idx] +=
                                    chain_data.get_edge_length(chain_eid) / edge_len_sum;
                            } else {
                                internal_costs[idx] +=
                                    chain_data.get_edge_length(chain_eid) / edge_len_sum;
                            }
                        }
                    }
                }
            }
        }

        let mut arc_builder = ArcBuilder::new(internal_costs, external_costs);
        for (eid, is_boundary) in edge_is_boundary_map {
            if is_boundary {
                continue;
            }

            let mut two_side_patch_node = [INVALID_IND; 2];
            for he in patch_mesh.edge(eid).halfedges() {
                if let Some(&idx) = patch_to_idx_map.get(&*he.face()) {
                    if he.same_dir() {
                        two_side_patch_node[0] = idx;
                    } else {
                        two_side_patch_node[1] = idx;
                    }
                    if two_side_patch_node.iter().all(|&x| x != INVALID_IND) {
                        break;
                    }
                }
            }

            if two_side_patch_node.iter().all(|&x| x != INVALID_IND) {
                arc_builder.add_arc(
                    two_side_patch_node[0],
                    two_side_patch_node[1],
                    chain_data.get_edge_length(eid) / edge_len_sum,
                    true,
                );
            }
        }

        let mut max_flow = PushRelabelFifo::from((arc_builder.arcs, patch_to_idx_map.len() + 2));
        max_flow.find_max_flow();

        Vec::from_iter(patch_to_idx_map.into_iter().filter_map(|(fid, idx)| {
            if max_flow.is_sink(idx + 1) {
                Some(fid)
            } else {
                None
            }
        }))
    }

    fn find_model_occupied_cells(
        &self,
        first_cid: usize,
        bdy_patch_masks: &mut [u8],
        cell_visited: &mut [bool],
    ) -> Option<Vec<usize>> {
        let mut occupied_cells = vec![first_cid];
        cell_visited[first_cid] = true;
        let mut idx = 0;
        let mut valid = true;
        loop {
            if idx >= occupied_cells.len() {
                break;
            }
            let curr_cid = occupied_cells[idx];
            idx += 1;
            for &oriented_pid in &self.cells[curr_cid] {
                let mask = bdy_patch_masks[strip_orientation(oriented_pid)];
                if mask == 2 {
                    let twin_cid = self.patch_cell_arr[twin_index(oriented_pid)];
                    if !cell_visited[twin_cid] {
                        cell_visited[twin_cid] = true;
                        occupied_cells.push(twin_cid);
                    }
                } else if (mask == 1) == is_positive(oriented_pid) {
                    valid = false;
                    break;
                }
            }
            if !valid {
                break;
            }
        }
        for &cid in &occupied_cells {
            cell_visited[cid] = false;
        }

        if valid { Some(occupied_cells) } else { None }
    }

    fn find_model_occupied_cells_fallback<T: IntoIterator<Item = usize>>(
        &self,
        model_bry_ori_patches: T,
    ) -> Vec<usize> {
        let mut internal_cost = vec![0.0; self.cells.len() + 1];
        let mut external_cost = vec![0.0; self.cells.len() + 1];
        let mut boundary_patches = HashSet::new();
        for ori_pid in model_bry_ori_patches {
            let pid = strip_orientation(ori_pid);
            boundary_patches.insert(pid);
            let area = self.get_patch_area(pid);
            let twin_ori_pid = twin_index(ori_pid);
            external_cost[self.patch_cell_arr[twin_ori_pid]] += area;
            internal_cost[self.patch_cell_arr[ori_pid]] += area;
        }
        internal_cost[self.cells.len()] = 1.0;

        let mut builder = ArcBuilder::new(internal_cost, external_cost);
        for pid in 0..self.patches.len() {
            if boundary_patches.contains(&pid) {
                continue;
            }
            let pos_ori_pid = oriented_index(pid, false);
            let neg_ori_pid = twin_index(pos_ori_pid);
            builder.add_arc(
                self.patch_cell_arr[pos_ori_pid],
                self.patch_cell_arr[neg_ori_pid],
                self.get_patch_area(pid),
                true,
            );
        }
        let mut max_flow = PushRelabelFifo::from((builder.arcs, self.cells.len() + 3));
        max_flow.find_max_flow();

        Vec::from_iter((0..self.cells.len()).filter_map(|cid| {
            if max_flow.is_sink(cid + 1) {
                Some(cid)
            } else {
                None
            }
        }))
    }

    fn build_patch_mesh(
        &self,
        chains: &[Vec<HalfedgeId>],
        surfaces: &[Surf],
        non_manifold_vertices: &[VertexId],
    ) -> (HoleAwareMesh<std::alloc::Global>, Vec<Bitmask>, Vec<usize>) {
        let mut patch_oriented_chains_arr = vec![Vec::with_capacity(4); self.patches.len()];
        for (chain_id, chain) in chains.iter().enumerate() {
            let chain_hid = chain[0];
            for he in self.mesh.halfedge(chain_hid).edge().halfedges() {
                let patch_id = self.face_patch_arr[*he.face()];
                patch_oriented_chains_arr[patch_id].push(oriented_index(
                    chain_id,
                    !self.mesh.hes_same_dir(chain_hid, *he),
                ));
            }
        }

        let non_manifold_vert_idx_map = HashMap::<_, _>::from_iter(
            non_manifold_vertices
                .iter()
                .enumerate()
                .map(|(idx, &vid)| (vid, idx)),
        );

        let mut chain_visited = vec![false; chains.len()];
        let mut loops = TwoDimArr::new();
        let mut face_loops = TwoDimArr::new();
        for (pid, ori_chains) in patch_oriented_chains_arr.into_iter().enumerate() {
            let surf = &surfaces[self.patch_surface_arr[pid]];
            let start = loops.len();
            for wire in self.get_patch_wires(chains, &ori_chains, &mut chain_visited, surf) {
                loops.push(wire.into_iter().map(|ori_chain_id| {
                    *non_manifold_vert_idx_map
                        .get(&get_ori_chain_vb(chains, ori_chain_id, &self.mesh))
                        .unwrap()
                }));
            }
            face_loops.push(start..loops.len());
        }
        let mesh = HoleAwareMesh::new(loops.iter(), face_loops.iter(), std::alloc::Global);
        let chain_masks = Vec::from_iter(mesh.edges().map(|edge| {
            let mut edge_mask = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
            for he in edge.halfedges() {
                edge_mask.set(self.patch_surface_arr[*he.face()]); //`*he.face()` means patch id
            }
            edge_mask
        }));
        let mut edge_chain_indices = vec![INVALID_IND; mesh.n_edges()];
        for (chain_id, chain) in chains.iter().enumerate() {
            let [va, vb] = get_chain_vertices(&chain, &self.mesh)
                .map(|vid| *non_manifold_vert_idx_map.get(&vid).unwrap());
            let eid = mesh.e_from_va_vb(va.into(), vb.into());
            edge_chain_indices[eid] = chain_id;
        }
        (mesh, chain_masks, edge_chain_indices)
    }

    fn get_patch_wires(
        &self,
        chains: &[Vec<HalfedgeId>],
        patch_oriented_chains: &[usize],
        chain_visited: &mut Vec<bool>,
        surf: &Surf,
    ) -> Vec<Vec<usize>> {
        let mut vertex_to_ori_chain_map = HashMap::<VertexId, TinyVec<[usize; 2]>>::new();
        for &oriented_chain_id in patch_oriented_chains {
            for vid in get_chain_vertices(&chains[strip_orientation(oriented_chain_id)], &self.mesh)
            {
                vertex_to_ori_chain_map
                    .entry(vid)
                    .or_default()
                    .push(oriented_chain_id);
            }
        }

        let mut wires = Vec::new();

        let propagate_wire = |first_ori_chain: usize, chain_visited: &mut [bool]| {
            let mut wire = vec![first_ori_chain];
            let mut current_ori_chain = first_ori_chain;
            let mut current_vid = get_ori_chain_vb(chains, current_ori_chain, &self.mesh);
            loop {
                let mut next_ori_chain = INVALID_IND;
                for &ori_chain_id in vertex_to_ori_chain_map.get(&current_vid).unwrap() {
                    if chain_visited[strip_orientation(ori_chain_id)] {
                        continue;
                    }
                    next_ori_chain = ori_chain_id;
                }

                if next_ori_chain == INVALID_IND {
                    break;
                }

                wire.push(next_ori_chain);
                current_ori_chain = next_ori_chain;
                chain_visited[strip_orientation(current_ori_chain)] = true;
                current_vid = get_ori_chain_vb(chains, current_ori_chain, &self.mesh);
            }
            wire
        };
        for &oriented_chain_id in patch_oriented_chains {
            let chain_id = strip_orientation(oriented_chain_id);
            if chain_visited[chain_id] {
                continue;
            }
            chain_visited[chain_id] = true;
            wires.push(propagate_wire(oriented_chain_id, chain_visited));
        }

        for &oriented_chain_id in patch_oriented_chains {
            chain_visited[strip_orientation(oriented_chain_id)] = false;
        }
        if wires.len() > 1 {
            let idx = wires
                .iter()
                .position(|wire| self.get_wire_uv_area(wire, chains, surf) > 0.0);
            if let Some(idx) = idx
                && idx > 0
            {
                wires.swap(0, idx);
            }
        }
        wires
    }

    fn get_wire_uv_area(&self, wire: &[usize], chains: &[Vec<HalfedgeId>], surf: &Surf) -> f64 {
        let mut uv_points = Vec::new();

        for &ori_chain_id in wire {
            if is_positive(ori_chain_id) {
                self.push_uv_stream(
                    &chains[strip_orientation(ori_chain_id)],
                    false,
                    surf,
                    &mut uv_points,
                );
            } else {
                self.push_uv_stream(
                    chains[strip_orientation(ori_chain_id)].iter().rev(),
                    true,
                    surf,
                    &mut uv_points,
                );
            }
        }
        face_area_2d(&uv_points)
    }

    fn push_uv_stream<'a, S: IntoIterator<Item = &'a HalfedgeId>>(
        &self,
        chain_stream: S,
        reversed: bool,
        surf: &Surf,
        uv_points: &mut Vec<f64>,
    ) {
        for &hid in chain_stream.into_iter() {
            let vid = if reversed {
                self.mesh.he_from(hid)
            } else {
                self.mesh.he_to(hid)
            };
            let pt = self.v_point(vid);
            uv_points.extend_from_slice(&if uv_points.is_empty() {
                surf.uv(pt, None)
            } else {
                surf.uv(pt, Some(&uv_points[(uv_points.len() - 2)..]))
            });
        }
    }

    fn propagate_edge_chain(
        &self,
        edge_mask: &Bitmask,
        va: VertexId,
        vb: VertexId,
        patch_mesh: &HoleAwareMesh<std::alloc::Global>,
        chain_masks: &[Bitmask],
        result: &mut Vec<HalfedgeId>,
    ) -> bool {
        if va == vb {
            return true;
        }
        for he in patch_mesh.vertex(va).outgoing_halfedges() {
            let eid = *he.edge();
            let mask = &chain_masks[eid];
            if let Some(&prev_hid) = result.last()
                && patch_mesh.he_edge(prev_hid) == eid
            {
                continue;
            }
            if mask.contain(edge_mask) {
                result.push(*he);
                if self.propagate_edge_chain(
                    edge_mask,
                    *he.to(),
                    vb,
                    patch_mesh,
                    chain_masks,
                    result,
                ) {
                    return true;
                } else {
                    result.pop();
                }
            }
        }
        false
    }

    #[inline]
    fn v_point(&self, vid: VertexId) -> &[f64] {
        point_3(&self.points, vid.0)
    }

    fn get_patch_area(&self, pid: usize) -> f64 {
        if let Some(patch_areas) = self.patch_areas.borrow().as_ref() {
            return patch_areas[pid];
        }
        let mut patch_areas = self
            .patches
            .iter()
            .map(|faces| {
                let mut area = 0.0;
                for &fid in faces {
                    let polygon = self
                        .mesh
                        .face(fid)
                        .vertices()
                        .map(|v| self.v_point(*v))
                        .collect_vec();
                    let pa = &polygon[0];
                    area += polygon[1..]
                        .windows(2)
                        .map(|pts| {
                            let v1 = sub_short::<3, _>(pts[0], pa);
                            let v2 = sub_short::<3, _>(pts[1], pa);
                            norm(&cross(&v1, &v2))
                        })
                        .sum::<f64>();
                }
                area
            })
            .collect_vec();
        let area_sum: f64 = patch_areas.iter().sum();
        patch_areas.iter_mut().for_each(|area| *area /= area_sum);
        let ret = patch_areas[pid];
        self.patch_areas.borrow_mut().replace(patch_areas);
        ret
    }

    fn get_result_oriented_patches(
        &self,
        model_cells_arr: Vec<Vec<usize>>,
        func: impl Fn(&[bool]) -> bool,
    ) -> Vec<usize> {
        let mut cell_is_in_model_interior =
            vec![vec![false; model_cells_arr.len()]; self.cells.len()];
        for (i, model_cells) in model_cells_arr.into_iter().enumerate() {
            for cid in model_cells {
                cell_is_in_model_interior[cid][i] = true;
            }
        }

        let mut is_cell_kept = Vec::with_capacity(cell_is_in_model_interior.len() + 1);

        is_cell_kept.extend(
            cell_is_in_model_interior
                .into_iter()
                .map(|flags| func(&flags)),
        );
        is_cell_kept.push(false);
        Vec::from_iter((0..self.patch_surface_arr.len()).filter_map(|pid| {
            let pos_pid = oriented_index(pid, false);
            let neg_pid = oriented_index(pid, true);
            let pos_cid = self.patch_cell_arr[pos_pid];
            let neg_cid = self.patch_cell_arr[neg_pid];
            if is_cell_kept[pos_cid] != is_cell_kept[neg_cid] {
                if is_cell_kept[pos_cid] {
                    Some(neg_pid)
                } else {
                    Some(pos_pid)
                }
            } else {
                None
            }
        }))
    }
}

struct ChainData<'a> {
    chains: Vec<Vec<HalfedgeId>>,
    patch_edge_chain_indices: Vec<usize>,
    chain_lengths: Vec<f64>,
    mesh: &'a SurfaceMesh,
    points: &'a [f64],
}

impl<'a> ChainData<'a> {
    fn new(
        chains: Vec<Vec<HalfedgeId>>,
        patch_edge_chain_indices: Vec<usize>,
        mesh: &'a SurfaceMesh,
        points: &'a [f64],
    ) -> Self {
        let chain_lengths = vec![f64::NAN; chains.len()];
        Self {
            chains,
            patch_edge_chain_indices,
            mesh,
            points,
            chain_lengths,
        }
    }

    fn get_edge_length(&mut self, patch_eid: EdgeId) -> f64 {
        let chain_id = self.patch_edge_chain_indices[patch_eid];
        if !self.chain_lengths[chain_id].is_nan() {
            return self.chain_lengths[chain_id];
        }
        let chain = &self.chains[chain_id];
        let mut prev_pt = point::<3>(self.points, *self.mesh.he_from(chain[0]));
        chain
            .iter()
            .map(|&hid| {
                let pt = point::<3>(self.points, *self.mesh.he_to(hid));
                let length = norm(&sub_short::<3, _>(pt, prev_pt));
                prev_pt = pt;
                length
            })
            .sum()
    }
}

fn get_chain_vertices<M: Mesh>(chain: &[HalfedgeId], mesh: &M) -> [VertexId; 2] {
    [get_chain_va(chain, mesh), get_chain_vb(chain, mesh)]
}

fn get_chain_va<M: Mesh>(chain: &[HalfedgeId], mesh: &M) -> VertexId {
    mesh.he_from(chain[0])
}

fn get_chain_vb<M: Mesh>(chain: &[HalfedgeId], mesh: &M) -> VertexId {
    mesh.he_to(chain[chain.len() - 1])
}

fn get_ori_chain_vertex<M: Mesh>(
    chains: &[Vec<HalfedgeId>],
    ori_chain_id: usize,
    is_first: bool,
    mesh: &M,
) -> VertexId {
    let chain_id = strip_orientation(ori_chain_id);
    let chain = &chains[chain_id];
    if is_positive(ori_chain_id) == is_first {
        get_chain_va(chain, mesh)
    } else {
        get_chain_vb(chain, mesh)
    }
}

fn get_ori_chain_vb<M: Mesh>(
    chains: &[Vec<HalfedgeId>],
    ori_chain_id: usize,
    mesh: &M,
) -> VertexId {
    get_ori_chain_vertex(chains, ori_chain_id, false, mesh)
}
