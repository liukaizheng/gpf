use std::{alloc::Allocator, collections::VecDeque};

use crate::{
    mesh::{ElementId, FaceId, HalfedgeId, ManifoldMesh, Mesh, VertexId},
    point_3,
    predicates::{Orientation, max_comp_in_tri_normal, miss_alignment, orient3d::orient3d_eeee},
};

use super::triangulate_points;

pub enum Convex3Result<A: Allocator + Copy> {
    Dim0(usize),
    Dim1(Vec<usize, A>),
    Dim2(Vec<usize, A>),
    Dim3(Vec<usize, A>),
}

pub fn convex_3<A: Allocator + Copy>(
    points: &[f64],
    only_verts: bool,
    alloc: A,
) -> Convex3Result<A> {
    let n_points = points.len() / 3;
    let mut indices = Vec::with_capacity_in(n_points, alloc);
    indices.extend(0..n_points);

    indices.sort_unstable_by(|&i, &j| {
        let pa = point_3(points, i);
        let pb = point_3(points, j);
        pa.partial_cmp(pb).unwrap()
    });

    indices.dedup_by(|i, j| {
        let pa = point_3(points, *i);
        let pb = point_3(points, *j);
        pa[0] == pb[0] && pa[1] == pb[1] && pa[2] == pb[2]
    });

    if indices.len() == 1 {
        return Convex3Result::Dim0(indices[0]);
    }

    let collinear_len = hull_1(points, &indices, alloc);
    if collinear_len == indices.len() {
        return Convex3Result::Dim1(indices);
    }

    let coplanar_len = hull_2(points, &indices, collinear_len, alloc);

    let triangles = if coplanar_len == 3 {
        let mut triangles = Vec::with_capacity_in(3, alloc);
        triangles.extend(&indices[..3]);
        triangles
    } else {
        let axis = max_comp_in_tri_normal(
            point_3(points, indices[0]),
            point_3(points, indices[1]),
            point_3(points, indices[collinear_len]),
            alloc,
        );

        let mut points_2d = Vec::new_in(alloc);

        points_2d.extend(indices[..coplanar_len].into_iter().flat_map(|&idx| {
            let p = point_3(points, idx);
            [p[(axis + 1) % 3], p[(axis + 2) % 3]]
        }));
        let mut triangles = Vec::new_in(alloc);
        triangles.extend(
            triangulate_points(&points_2d, true, alloc)
                .into_iter()
                .map(|idx| indices[idx]),
        );

        triangles
    };
    if coplanar_len == indices.len() {
        return Convex3Result::Dim2(triangles);
    }
    let mesh = hull_3(points, &indices, coplanar_len, triangles, alloc);
    let result = if only_verts {
        let mut result = Vec::with_capacity_in(mesh.n_vertices(), alloc);
        result.extend(mesh.vertices().map(|v| v.0));
        result
    } else {
        let mut result = Vec::with_capacity_in(mesh.n_faces(), alloc);
        result.extend(mesh.faces().flat_map(|f| {
            let he = f.halfedge();
            let he_next = he.next();
            let he_nnext = he_next.next();
            [he.to().0, he_next.to().0, he_nnext.to().0]
        }));
        result
    };

    Convex3Result::Dim3(result)
}

fn hull_1<A: Allocator + Copy>(points: &[f64], indices: &[usize], alloc: A) -> usize {
    let pa = point_3(points, indices[0]);
    let pb = point_3(points, indices[1]);
    for i in 2..indices.len() {
        let pc = point_3(points, indices[i]);
        if miss_alignment(pa, pb, pc, alloc) {
            return i;
        }
    }
    return indices.len();
}

fn hull_2<A: Allocator + Copy>(points: &[f64], indices: &[usize], start: usize, alloc: A) -> usize {
    let pa = point_3(points, indices[0]);
    let pb = point_3(points, indices[1]);
    let pc = point_3(points, indices[start]);
    for i in (start + 1)..indices.len() {
        let pd = point_3(points, indices[i]);
        if orient3d_eeee(pa, pb, pc, pd, alloc) != Orientation::Zero {
            return i;
        }
    }
    return indices.len();
}

fn hull_3<A: Allocator + Copy>(
    points: &[f64],
    indices: &[usize],
    start: usize,
    mut triangles: Vec<usize, A>,
    alloc: A,
) -> ManifoldMesh<A> {
    let is_neg = {
        let pa = point_3(points, triangles[0]);
        let pb = point_3(points, triangles[1]);
        let pc = point_3(points, triangles[2]);
        let pd = point_3(points, indices[start]);
        orient3d_eeee(pa, pb, pc, pd, alloc) == Orientation::Positive
    };

    if is_neg {
        triangles.reverse();
    }
    let mut mesh = ManifoldMesh::new(triangles.chunks(3), alloc);
    mesh.set_n_vertices(points.len() / 3);

    let mut first_bot_hid = HalfedgeId::default();
    for v in mesh.vertices() {
        let hid = *v.out_halfedge();
        if mesh.he_is_boundary(hid) {
            first_bot_hid = hid;
            break;
        }
    }
    close_hull(&mut mesh, first_bot_hid, indices[start].into());

    let mut visited = Vec::new_in(alloc);
    visited.resize(mesh.n_faces_capacity(), false);
    let mut kept = Vec::new_in(alloc);
    kept.resize(mesh.n_faces_capacity(), false);

    let mut visible_faces = VecDeque::new_in(alloc);
    let mut removed_faces = Vec::new_in(alloc);
    let mut kept_faces = Vec::new_in(alloc);
    let mut prev_vid = indices[start].into();
    for &pid in &indices[(start + 1)..] {
        let pd = point_3(points, pid);
        let vid = pid.into();
        let visible = |fid: FaceId| {
            let he = mesh.face(fid).halfedge();

            let pa = point_3(points, he.to().0);
            let he = he.next();
            let pb = point_3(points, he.to().0);
            let he = he.next();
            let pc = point_3(points, he.to().0);
            orient3d_eeee(pa, pb, pc, pd, alloc) == Orientation::Positive
        };
        visible_faces.clear();
        removed_faces.clear();
        kept_faces.clear();
        visible_faces.extend(mesh.vertex(prev_vid).incoming_halfedges().filter_map(|he| {
            let fid = mesh.he_face(*he);
            debug_assert!(fid.valid());
            visited[fid.0] = true;
            if visible(fid) {
                removed_faces.push(fid);
                Some(fid)
            } else {
                kept[fid.0] = true;
                kept_faces.push(fid);
                None
            }
        }));
        debug_assert!(!visible_faces.is_empty());

        first_bot_hid = HalfedgeId::default();
        loop {
            if visible_faces.is_empty() {
                break;
            }
            let fid = visible_faces.pop_front().unwrap();
            for he in mesh.face(fid).halfedges() {
                let twin_he = he.twin();
                let twin_fid = mesh.he_face(*twin_he);
                if visited[twin_fid.0] {
                    continue;
                }

                visited[twin_fid.0] = true;

                if visible(twin_fid) {
                    visible_faces.push_back(twin_fid);
                    removed_faces.push(twin_fid);
                } else {
                    kept_faces.push(twin_fid);
                    kept[twin_fid.0] = true;
                    if !first_bot_hid.valid() {
                        first_bot_hid = *he;
                    }
                }
            }
        }
        if !first_bot_hid.valid() {
            for &fid in &kept_faces {
                for he in mesh.face(fid).halfedges() {
                    let twin_hid = *he.twin();
                    let twin_fid = mesh.he_face(twin_hid);
                    if visited[twin_fid.0] && !kept[twin_fid.0] {
                        first_bot_hid = twin_hid;
                        break;
                    }
                }
                if first_bot_hid.valid() {
                    break;
                }
            }
        }

        debug_assert!(first_bot_hid.valid());
        for &fid in &removed_faces {
            mesh.remove_face(fid);
        }
        debug_assert!(mesh.he_is_boundary(first_bot_hid));

        for &fid in &kept_faces {
            visited[fid.0] = false;
            kept[fid.0] = false;
        }

        close_hull(&mut mesh, first_bot_hid, vid);

        visited.resize(mesh.n_faces_capacity(), false);
        kept.resize(mesh.n_faces_capacity(), false);
        prev_vid = vid;
    }
    mesh
}

fn close_hull<A: Allocator + Copy>(
    mesh: &mut ManifoldMesh<A>,
    first_hid: HalfedgeId,
    vid: VertexId,
) {
    debug_assert!(mesh.he_is_boundary(first_hid));
    let start_vid = mesh.he_from(first_hid);

    let mut curr_bot_hid = first_hid;
    let mut first_side_hid = HalfedgeId::default();
    let mut prev_side_hid = HalfedgeId::default();
    loop {
        let mut next_bot_hid = mesh.he_next(curr_bot_hid);
        if mesh.he_to(curr_bot_hid) == start_vid {
            next_bot_hid = first_hid;
        }
        let [va, vb] = mesh.he_vertices(curr_bot_hid);
        if curr_bot_hid == first_hid {
            prev_side_hid = mesh.new_edge_by_veritces(vid, va);
            first_side_hid = mesh.he_twin(prev_side_hid);
        }

        let curr_side_hid = if next_bot_hid == first_hid {
            first_side_hid
        } else {
            mesh.new_edge_by_veritces(vb, vid)
        };

        mesh.new_face_by_halfedges(&[prev_side_hid, curr_bot_hid, curr_side_hid]);
        if next_bot_hid == first_hid {
            break;
        }
        curr_bot_hid = next_bot_hid;
        prev_side_hid = mesh.he_twin(curr_side_hid);
    }
}
