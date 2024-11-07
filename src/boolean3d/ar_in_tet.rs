use std::alloc::{Allocator, Global};

use bumpalo::Bump;
use std::collections::HashMap;
use tinyvec::{tiny_vec, TinyVec};

use crate::{
    math::interpolate,
    mesh::{clone_vec_in, EdgeId, ElementId, FaceId, Mesh, SurfaceMesh, VertexId},
    point,
    predicates::det4,
    signed_index, INVALID_IND,
};

use super::tet_set::{tet_face_reversed, TetSet};

#[derive(Clone)]
enum InterPt {
    V(VertexId),
    ES((EdgeId, usize)),
    FSS((FaceId, usize, usize)),
    SSS([usize; 3]),
    INVALID,
}

impl InterPt {
    #[inline]
    fn valid(&self) -> bool {
        match self {
            InterPt::INVALID => false,
            _ => true,
        }
    }
}

#[derive(Clone)]
struct VertexData {
    planes: [usize; 3],
    parents: [VertexId; 2],
    vals: [f64; 2],
}

#[derive(Clone)]
struct FaceData {
    /// plane id
    pid: usize,
    /// `cells[0]`: inner cell of this face
    /// `cells[1]`: outer cell of this face
    cells: [usize; 2],

    surfaces: TinyVec<[i64; 1]>,
}

pub(super) struct Arrangement<A: Allocator + Copy> {
    mesh: SurfaceMesh<A>,
    vertices: Vec<VertexData, A>,
    edges: Vec<[usize; 2], A>,
    planes: Vec<[f64; 4], A>,
    face_data: Vec<FaceData, A>,
    cell_faces: Vec<Vec<FaceId, A>, A>,
}

fn orient3d<A: Allocator + Copy>(
    pa: &[f64; 4],
    pb: &[f64; 4],
    pc: &[f64; 4],
    pd: &[f64; 4],
    alloc: A,
) -> f64 {
    #[rustfmt::skip]
    let numerator = det4(
        pa[0], pa[1], pa[2], pa[3],
        pb[0], pb[1], pb[2], pb[3],
        pc[0], pc[1], pc[2], pc[3],
        pd[0], pd[1], pd[2], pd[3],
        alloc,
    );
    #[rustfmt::skip]
    let denominator = det4(
        pa[0], pa[1], pa[2], pa[3],
        pb[0], pb[1], pb[2], pb[3],
        pc[0], pc[1], pc[2], pc[3],
        1.0  , 1.0  , 1.0  , 1.0  ,
        alloc,
    );

    debug_assert!(denominator != 0.0);
    if denominator > 0.0 {
        return numerator;
    } else {
        return -numerator;
    }
}

impl<A: Allocator + Copy> Arrangement<A> {
    fn new_tet(alloc: A) -> Self {
        let mesh = SurfaceMesh::new([[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]], alloc);
        let mut vertices = Vec::with_capacity_in(4, alloc);
        vertices.extend(mesh.vertices().map(|v| {
            let mut planes = [0; 3];
            for (fid, hid) in planes.iter_mut().zip(v.incoming_halfedges()) {
                *fid = mesh.he_face(*hid).0;
            }
            VertexData {
                planes,
                parents: [VertexId::default(); 2],
                vals: [f64::NAN; 2],
            }
        }));

        let mut edges = Vec::with_capacity_in(6, alloc);
        edges.extend(mesh.edges().map(|e| {
            let mut res = [0; 2];
            for (fid, hid) in res.iter_mut().zip(e.halfedges()) {
                *fid = mesh.he_face(*hid).0;
            }
            res
        }));

        let mut planes = Vec::with_capacity_in(4, alloc);
        planes.push([1.0, 0.0, 0.0, 0.0]);
        planes.push([0.0, 1.0, 0.0, 0.0]);
        planes.push([0.0, 0.0, 1.0, 0.0]);
        planes.push([0.0, 0.0, 0.0, 1.0]);

        let mut face_data = Vec::with_capacity_in(4, alloc);
        face_data.extend((0..4).map(|pid| FaceData {
            pid,
            cells: [0, INVALID_IND],
            surfaces: TinyVec::new(),
        }));

        let mut one_cell_faces = Vec::with_capacity_in(4, alloc);
        one_cell_faces.extend((0..4).map(|idx| FaceId::from(idx)));

        let mut cell_faces = Vec::with_capacity_in(1, alloc);
        cell_faces.push(one_cell_faces);

        Self {
            mesh,
            vertices,
            edges,
            planes,
            face_data,
            cell_faces,
        }
    }

    fn clone_in<A1: Allocator + Copy>(&self, alloc: A1) -> Arrangement<A1> {
        let mut cell_faces = Vec::with_capacity_in(self.cell_faces.len(), alloc);
        for v in &self.cell_faces {
            cell_faces.push(clone_vec_in(v, alloc));
        }
        Arrangement {
            mesh: self.mesh.clone_in(alloc),
            vertices: clone_vec_in(&self.vertices, alloc),
            edges: clone_vec_in(&self.edges, alloc),
            planes: clone_vec_in(&self.planes, alloc),
            face_data: clone_vec_in(&self.face_data, alloc),
            cell_faces,
        }
    }

    #[inline]
    fn is_face_inner_cell(&self, fid: FaceId, cid: usize) -> bool {
        let cells = &self.face_data[fid].cells;
        debug_assert!(cells[0] == cid || cells[1] == cid);
        cells[0] == cid
    }

    fn add_plane<A1: Allocator + Copy>(&mut self, pid: usize, sid: usize, alloc: A1) {
        let mut vert_orientations = Vec::with_capacity_in(self.mesh.n_vertices_capacity(), alloc);
        vert_orientations.extend(self.mesh.vertices().map(|v| {
            let vid = *v;
            let v_planes = &self.vertices[vid].planes;
            orient3d(
                &self.planes[v_planes[0]],
                &self.planes[v_planes[1]],
                &self.planes[v_planes[2]],
                &self.planes[pid],
                alloc,
            )
        }));

        let mut coplanar_pid = INVALID_IND;
        let mut non_zero_ori = f64::NAN;
        for face in self.mesh.faces() {
            let fid = *face;
            let face_pid = self.face_data[fid].pid;
            if face_pid < 4 {
                continue;
            }

            if coplanar_pid != INVALID_IND {
                debug_assert!(!non_zero_ori.is_nan());
                if face_pid == coplanar_pid {
                    self.face_data[fid]
                        .surfaces
                        .push(signed_index(sid, non_zero_ori > 0.0));
                }
            } else {
                if face.vertices().all(|v| vert_orientations[*v] == 0.0) {
                    // println!("add same plane");
                    let base_fid = *face;
                    let pid = self.face_data[base_fid].pid;
                    debug_assert!(pid >= 4);
                    let pos_cid = self.face_data[pid].cells[0];
                    debug_assert!(pos_cid != INVALID_IND);

                    let non_zero_ori_fn = || {
                        for &fid in &self.cell_faces[pos_cid] {
                            if fid == base_fid {
                                continue;
                            }

                            for v in self.mesh.face(fid).vertices() {
                                let ori = vert_orientations[*v];
                                if ori != 0.0 {
                                    return ori;
                                }
                            }
                        }
                        return f64::NAN;
                    };
                    coplanar_pid = face_pid;
                    non_zero_ori = non_zero_ori_fn();
                    debug_assert!(!non_zero_ori.is_nan());
                    self.face_data[fid]
                        .surfaces
                        .push(signed_index(sid, non_zero_ori > 0.0));
                }
            }
        }

        if coplanar_pid != INVALID_IND {
            return;
        }

        self.split_edges(&mut vert_orientations, pid, alloc);
        self.split_faces(&vert_orientations, pid);
        self.split_cells(&vert_orientations, pid, sid, alloc);
    }

    fn split_edges<A1: Allocator + Copy>(
        &mut self,
        vert_orientations: &mut Vec<f64, A1>,
        pid: usize,
        alloc: A1,
    ) {
        let n_old_edges = self.mesh.n_edges_capacity();
        for eid in 0..n_old_edges {
            let eid = eid.into();
            let [va, vb] = self.mesh.e_vertices(eid);
            let ori1 = vert_orientations[va];
            let ori2 = vert_orientations[vb];
            if (ori1 < 0.0 && ori2 > 0.0) || (ori1 > 0.0 && ori2 < 0.0) {
                self.mesh.split_edge(eid, alloc);
                let e_planes = &self.edges[eid];
                self.vertices.push(VertexData {
                    planes: [e_planes[0], e_planes[1], pid],
                    parents: [va, vb],
                    vals: [ori1, ori2],
                });
                vert_orientations.push(0.0);
                self.edges.push(e_planes.clone());
            }
        }
    }

    fn split_faces(&mut self, vert_orientations: &[f64], pid: usize) {
        let n_old_faces = self.mesh.n_faces_capacity();
        for fid in 0..n_old_faces {
            let fid = fid.into();
            let mut first_zero_vid = VertexId::default();
            let first_hid = self.mesh.f_halfedge(fid);
            let mut curr_hid = first_hid;
            let mut prev_ori = vert_orientations[self.mesh.he_from(curr_hid)];
            if prev_ori == 0.0 {
                first_zero_vid = self.mesh.he_from(curr_hid);
            }
            loop {
                let vid = self.mesh.he_to(curr_hid);
                let ori = vert_orientations[vid];
                let next_hid = self.mesh.he_next(curr_hid);
                if ori == 0.0 {
                    if prev_ori == 0.0 {
                        self.mesh.set_f_halfedge(fid, curr_hid);
                        break;
                    } else {
                        if first_zero_vid.valid() {
                            if self.mesh.he_to(next_hid) == first_zero_vid {
                                self.mesh.set_f_halfedge(fid, next_hid);
                                break;
                            } else if vid == first_zero_vid {
                                break;
                            }
                            self.mesh.split_face(fid, first_zero_vid, vid);
                            let new_fid = self.face_data.len().into();
                            let data = self.face_data[fid].clone();
                            self.edges.push([pid, data.pid]);
                            for cid in data.cells {
                                if cid != INVALID_IND {
                                    self.cell_faces[cid].push(new_fid);
                                }
                            }
                            self.face_data.push(data);
                            break;
                        } else {
                            first_zero_vid = vid;
                        }
                    }
                }
                curr_hid = next_hid;
                if curr_hid == first_hid {
                    break;
                }

                prev_ori = ori;
            }
        }
    }

    fn split_cells<A1: Allocator + Copy>(
        &mut self,
        vert_orientations: &[f64],
        pid: usize,
        sid: usize,
        alloc: A1,
    ) {
        let n_old_cells = self.cell_faces.len();
        let self_alloc = self.vertices.allocator().clone();
        for cid in 0..n_old_cells {
            let mut pos_cell_faces = Vec::new_in(self_alloc);
            let mut neg_cell_faces = Vec::new_in(self_alloc);

            let mut start_zero_vid = VertexId::default();
            let mut n_halfedges = 0;
            for &fid in &self.cell_faces[cid] {
                let hid = self.mesh.f_halfedge(fid);
                let [va, vb] = self.mesh.he_vertices(hid);
                let mut non_zero_vid = VertexId::default();
                if vert_orientations[va] == 0.0 {
                    if vert_orientations[vb] == 0.0 {
                        let vc = *self.mesh.halfedge(hid).next().to();
                        debug_assert!(vert_orientations[vc] != 0.0);

                        let is_pos = vert_orientations[vc] > 0.0;
                        if is_pos {
                            pos_cell_faces.push(fid);
                        } else {
                            neg_cell_faces.push(fid);
                        }

                        if is_pos == self.is_face_inner_cell(fid, cid) {
                            self.mesh.set_v_halfedge(va, hid);
                            if !start_zero_vid.valid() {
                                start_zero_vid = va;
                            }
                            n_halfedges += 1;
                        }
                    } else {
                        non_zero_vid = vb;
                    }
                } else {
                    non_zero_vid = va;
                }

                if non_zero_vid.valid() {
                    if vert_orientations[non_zero_vid] > 0.0 {
                        pos_cell_faces.push(fid);
                    } else {
                        neg_cell_faces.push(fid);
                    }
                }
            }

            if n_halfedges < 3 {
                continue;
            }
            debug_assert!(start_zero_vid.valid());

            let mut curr_vid = start_zero_vid;
            let mut new_halfedges = Vec::with_capacity_in(n_halfedges, alloc);
            loop {
                let hid = self.mesh.v_halfedge(curr_vid);
                new_halfedges.push(hid);
                curr_vid = self.mesh.he_to(hid);
                if curr_vid == start_zero_vid {
                    break;
                }
            }

            let new_fid = self.mesh.add_face_by_halfedges(&new_halfedges);
            let new_cid = self.cell_faces.len();
            self.face_data.push(FaceData {
                pid,
                cells: [cid, new_cid],
                surfaces: tiny_vec!([i64; 1] => signed_index(sid, false)),
            });

            for &fid in &pos_cell_faces {
                for c in &mut self.face_data[fid].cells {
                    if *c == cid {
                        *c = new_cid;
                    }
                }
            }

            neg_cell_faces.push(new_fid);
            pos_cell_faces.push(new_fid);

            debug_assert!(neg_cell_faces.len() >= 4);
            debug_assert!(pos_cell_faces.len() >= 4);

            self.cell_faces[cid] = neg_cell_faces;
            self.cell_faces.push(pos_cell_faces);
        }
    }

    fn get_global_vertex(
        &self,
        vid: VertexId,
        tid: usize,
        tets: &TetSet,
        surfs: &[usize],
    ) -> InterPt {
        let mut boundries = TinyVec::<[usize; 3]>::new();
        let mut inners = TinyVec::<[usize; 3]>::new();
        for &i in &self.vertices[vid].planes {
            if i < 4 {
                boundries.push(i);
            } else {
                inners.push(i - 4);
            }
        }

        match inners.len() {
            0 => {
                let idx = boundries[0] ^ boundries[1] ^ boundries[2];
                InterPt::V(tets.tet_vertices[tid][idx])
            }
            1 => {
                let edge_index =
                    if boundries[0] != 0 { 0 } else { 1 } + 5 - boundries[0] - boundries[1];
                InterPt::ES((tets.tet_edges[tid][edge_index], surfs[inners[0]]))
            }
            2 => {
                let fid = tets.tet_faces[tid][boundries[0]];
                InterPt::FSS((fid, surfs[inners[0]], surfs[inners[1]]))
            }
            3 => InterPt::SSS([surfs[inners[0]], surfs[inners[1]], surfs[inners[2]]]),
            _ => InterPt::INVALID,
        }
    }

    fn extract_mesh(&mut self, tets: &TetSet, tid: usize, surfs: &[usize], data: &mut ExtractMesh) {
        let alloc = self.vertices.allocator();
        let mut vertex_pts = Vec::with_capacity_in(self.mesh.n_vertices_capacity(), alloc);
        vertex_pts.resize(self.mesh.n_vertices_capacity(), InterPt::INVALID);

        for face in self.mesh.faces() {
            let fid = *face;
            if self.face_data[fid].surfaces.is_empty() {
                continue;
            }
            for v in face.vertices() {
                let vid = *v;
                if !vertex_pts[vid].valid() {
                    vertex_pts[vid] = self.get_global_vertex(vid, tid, tets, surfs);
                }
            }
        }

        let mut vertex_indices = Vec::with_capacity_in(self.mesh.n_vertices_capacity(), alloc);
        vertex_indices.resize(self.mesh.n_vertices_capacity(), INVALID_IND);
        for (idx, inter_pt) in vertex_pts.into_iter().enumerate() {
            let pid = match inter_pt {
                InterPt::V(vid) => {
                    if data.point_map[vid] == INVALID_IND {
                        let pid = data.points.len() / 3;
                        data.points.extend_from_slice(point(&tets.points, vid.0));
                        pid
                    } else {
                        data.point_map[vid]
                    }
                }
                InterPt::ES(key) => match data.edge_point_map.entry(key) {
                    std::collections::hash_map::Entry::Occupied(occupied) => *occupied.get(),
                    std::collections::hash_map::Entry::Vacant(vacant) => {
                        let pid = data.points.len() / 3;
                        vacant.insert(pid);
                        pid
                    }
                },
                InterPt::FSS(key) => match data.face_point_map.entry(key) {
                    std::collections::hash_map::Entry::Occupied(occupied) => *occupied.get(),
                    std::collections::hash_map::Entry::Vacant(vacant) => {
                        let pid = data.points.len() / 3;
                        vacant.insert(pid);
                        pid
                    }
                },
                InterPt::SSS(_) => data.points.len() / 3,

                InterPt::INVALID => INVALID_IND,
            };

            if pid != INVALID_IND {
                if pid * 3 >= data.points.len() {
                    let [pa, pb] = self.vertices[idx].parents.map(|vid| {
                        if vid.0 < 4 {
                            point(&tets.points, tets.tet_vertices[tid][vid].0)
                        } else {
                            point(&data.points, vertex_indices[vid])
                        }
                    });
                    let [s1, s2] = self.vertices[idx].vals.map(|x| x.abs());
                    data.points
                        .extend_from_slice(&interpolate::<3>(pa, s1, pb, s2));
                }
                vertex_indices[idx] = pid;
            }
        }

        for face in self.mesh.faces() {
            let fid = *face;
            if self.face_data[fid].surfaces.is_empty() {
                continue;
            }

            let first_hid = self.mesh.f_halfedge(fid);
            let mut curr_hid = first_hid;
            let v1 = vertex_indices[self.mesh.he_to(curr_hid)];
            curr_hid = self.mesh.he_next(curr_hid);
            let mut v2 = vertex_indices[self.mesh.he_to(curr_hid)];
            loop {
                curr_hid = self.mesh.he_next(curr_hid);
                if curr_hid == first_hid {
                    break;
                }
                let v3 = vertex_indices[self.mesh.he_to(curr_hid)];
                data.triangles.extend_from_slice(&[v1, v2, v3]);
                data.triangle_parents
                    .push(self.face_data[fid].surfaces.clone());
                v2 = v3;
            }
        }
    }
}

struct ExtractMesh {
    points: Vec<f64>,
    triangles: Vec<usize>,
    triangle_parents: Vec<TinyVec<[i64; 1]>>,
    point_map: Vec<usize>,
    edge_point_map: HashMap<(EdgeId, usize), usize>,
    face_point_map: HashMap<(FaceId, usize, usize), usize>,
}

impl ExtractMesh {
    fn new(n_vertices: usize) -> Self {
        Self {
            points: Vec::new(),
            triangles: Vec::new(),
            triangle_parents: Vec::new(),
            point_map: vec![INVALID_IND; n_vertices],
            edge_point_map: HashMap::new(),
            face_point_map: HashMap::new(),
        }
    }
}

fn write_obj(name: &str, points: &[f64], triangles: &[usize]) {
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

    for i in 0..triangles.len() / 3 {
        writeln!(
            &mut file,
            "f {} {} {}",
            triangles[i * 3] + 1,
            triangles[i * 3 + 1] + 1,
            triangles[i * 3 + 2] + 1
        )
        .unwrap();
    }
}

pub(crate) fn extract_mesh(tets: &TetSet, vals: Vec<Vec<f64>>, active_surfaces: Vec<Vec<usize>>) {
    let base_tet = Arrangement::new_tet(Global);
    let mut data = ExtractMesh::new(tets.mesh.n_vertices_capacity());

    let mut bump = Bump::new();
    for (tid, verts) in tets.tet_vertices.iter().enumerate() {
        let surfs = &active_surfaces[tid];
        if surfs.is_empty() {
            continue;
        }
        bump.reset();

        let mut ar = base_tet.clone_in(&bump);
        ar.planes.reserve(ar.planes.len() + surfs.len());
        ar.planes.extend(surfs.iter().map(|&sid| {
            let mut tet_vals = [f64::NAN; 4];
            for (val, &vid) in tet_vals.iter_mut().zip(verts) {
                *val = vals[sid][vid];
            }
            tet_vals
        }));

        for (i, &fid) in tets.tet_faces[tid].iter().enumerate() {
            let face_surfs = &tets.face_surfaces[fid];
            if face_surfs.is_empty() {
                continue;
            }
            if tet_face_reversed(&tets.face_tets[fid], tid) {
                ar.face_data[i]
                    .surfaces
                    .extend(face_surfs.iter().map(|&sid| sid));
            } else {
                ar.face_data[i]
                    .surfaces
                    .extend(face_surfs.iter().map(|&sid| -sid));
            }
        }

        for (i, &sid) in surfs.iter().enumerate() {
            ar.add_plane(i + 4, sid, &bump);
        }

        ar.extract_mesh(tets, tid, surfs, &mut data);
    }

    write_obj("123.obj", &data.points, &data.triangles);
}
