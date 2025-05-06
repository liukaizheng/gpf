use hashbrown::HashMap;
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::geometry::{Surf, Surface};
use crate::math::{square_norm, sub_short};
use crate::mesh::{EdgeId, HalfedgeId};
use crate::{INVALID_IND, face_area_2d, is_positive, oriented_index, point_3, strip_orientation};
use crate::{
    boolean3d::{extract_cells::write_chains, write_obj},
    mesh::{FaceId, Mesh, SurfaceMesh, VertexId},
    utils::Bitmask,
};

use super::{BrepModel, extract_cells::identify_chain_edge};

pub(crate) struct ModelData {
    pub(crate) points: Vec<f64>,
    pub(crate) patches: Vec<Vec<FaceId>>,
    pub(crate) patch_surface_arr: Vec<usize>,
    pub(crate) surface_patches: Vec<Vec<usize>>,
    pub(crate) face_patch_arr: Vec<usize>,
    pub(crate) cells: Vec<Vec<usize>>,
    pub(crate) patch_cell_arr: Vec<usize>,
    pub(crate) mesh: SurfaceMesh,
}

impl ModelData {
    pub(crate) fn resolve(&self, models: Vec<BrepModel>, surfaces: &[Surf]) {
        write_obj("123.obj", &self.points, &self.mesh);
        let (non_manifold_vertices, chains, edge_chain_indices) =
            identify_chain_edge(&self.mesh, |fa, fb| {
                self.face_patch_arr[fa] == self.face_patch_arr[fb]
            });

        println!(
            "the number of non_manifold_vertices: {:?}",
            non_manifold_vertices.len()
        );
        println!("the number of chains: {:?}", chains.len());
        write_chains("chain.obj", &self.points, &self.mesh, &edge_chain_indices);

        let mut mask_vertices_map =
            HashMap::<Bitmask, TinyVec<[VertexId; 1]>>::with_capacity(non_manifold_vertices.len());

        let mut add_into_mask_vertices_map = |mask: Bitmask, vid: VertexId| {
            mask_vertices_map.entry(mask).or_default().push(vid);
        };

        for vid in non_manifold_vertices {
            let mut mask = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
            for f in self
                .mesh
                .vertex(vid)
                .incoming_halfedges()
                .map(|he| he.face())
            {
                mask.set(self.patch_surface_arr[self.face_patch_arr[*f]]);
            }

            if mask.n_elements() <= 3 {
                add_into_mask_vertices_map(mask, vid);
            } else {
                for k in 0..mask.n_elements() {
                    let mut m = Bitmask::<[usize; 1]>::new(self.surface_patches.len());
                    for elements in mask.iter_set_bits().combinations(k) {
                        m.set_from_iter(elements);
                    }
                    add_into_mask_vertices_map(m, vid);
                }
            }
        }

        let model_isomesh_vertices_map = models
            .iter()
            .map(|model| {
                model
                    .mesh
                    .vertices()
                    .map(|m_vert| {
                        let m_vid = *m_vert;
                        let mask = model.vert_mask(m_vid, self.surface_patches.len());
                        let mut iso_vid = VertexId::default();
                        let mut min_dist = f64::MAX;
                        let m_pt = model.v_point(m_vid);
                        for &vid in mask_vertices_map.get(&mask).unwrap_or(&TinyVec::new()) {
                            let pt = self.v_point(vid);
                            let dist = square_norm(&sub_short::<3, _>(m_pt, pt));
                            if dist < min_dist {
                                min_dist = dist;
                                iso_vid = vid;
                            }
                        }
                        iso_vid
                    })
                    .collect_vec()
            })
            .collect_vec();

        println!(
            "model_isomesh_vertices_map: {:?}",
            model_isomesh_vertices_map
        );
    }

    fn resolve_face_patches(
        &self,
        surfaces: &[Surf],
        model: &BrepModel,
        mask_vertices_map: &HashMap<Bitmask, TinyVec<[VertexId; 1]>>,
    ) {
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
                    let pt = self.v_point(vid);
                    let dist = square_norm(&sub_short::<3, _>(m_pt, pt));
                    if dist < min_dist {
                        min_dist = dist;
                        iso_vid = vid;
                    }
                }
                iso_vid
            })
            .collect_vec();
        for edge in model.mesh.edges() {
            self.propagate_edge_chain(*edge, model, &iso_vertices);
        }
    }

    fn build_patch_mesh(&self, chains: &[Vec<HalfedgeId>], surfaces: &[Surf]) {
        let mut patch_oriented_chain_arr = vec![Vec::with_capacity(4); self.patches.len()];
        for (chain_id, chain) in chains.iter().enumerate() {
            let chain_hid = chain[0];
            for he in self.mesh.halfedge(chain_hid).edge().halfedges() {
                let patch_id = self.face_patch_arr[*he.face()];
                patch_oriented_chain_arr[patch_id].push(oriented_index(
                    chain_id,
                    !self.mesh.hes_same_dir(chain_hid, *he),
                ));
            }
        }

        let mut chain_visited = vec![false; chains.len()];
        for (pid, ori_chains) in patch_oriented_chain_arr.into_iter().enumerate() {
            let surf = &surfaces[self.patch_surface_arr[pid]];
            let wires = self.get_patch_wires(chains, &ori_chains, &mut chain_visited, surf);
        }
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

        let propagate_wire = |first_ori_chain: usize, chain_visited: &mut Vec<bool>| {
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

    fn propagate_edge_chain(&self, eid: EdgeId, model: &BrepModel, iso_vertices: &[VertexId]) {
        let [va, vb] = model.mesh.e_vertices(eid).map(|vid| iso_vertices[vid]);
    }

    #[inline]
    fn v_point(&self, vid: VertexId) -> &[f64] {
        point_3(&self.points, vid.0)
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

fn get_ori_chain_va<M: Mesh>(
    chains: &[Vec<HalfedgeId>],
    ori_chain_id: usize,
    mesh: &M,
) -> VertexId {
    get_ori_chain_vertex(chains, ori_chain_id, true, mesh)
}

fn get_ori_chain_vb<M: Mesh>(
    chains: &[Vec<HalfedgeId>],
    ori_chain_id: usize,
    mesh: &M,
) -> VertexId {
    get_ori_chain_vertex(chains, ori_chain_id, false, mesh)
}
