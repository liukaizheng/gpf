use std::{alloc::Allocator, collections::VecDeque};

use crate::{
    mesh::{ElementId, FaceId, HalfedgeId, ManifoldMesh, Mesh, VertexId},
    point,
    predicates::{max_comp_in_tri_normal, orient2d, orient3d},
};

use super::triangulate;

pub enum Convex3Result<A: Allocator + Copy> {
    Dim0(usize),
    Dim1(Vec<usize, A>),
    Dim2(Vec<usize, A>),
    Dim3(Vec<usize, A>),
}

pub fn convex_3<A: Allocator + Copy>(points: &[f64], alloc: A) -> Convex3Result<A> {
    let n_points = points.len() / 3;
    let mut indices = Vec::with_capacity_in(n_points, alloc);
    indices.extend(0..n_points);

    indices.sort_unstable_by(|&i, &j| {
        let pa = point(points, i);
        let pb = point(points, j);
        pa.partial_cmp(pb).unwrap()
    });

    indices.dedup_by(|i, j| {
        let pa = point(points, *i);
        let pb = point(points, *j);
        pa[0] == pb[0] && pa[1] == pb[1] && pa[2] == pb[2]
    });

    if indices.len() == 1 {
        return Convex3Result::Dim0(indices[0]);
    }

    let colinear_len = hull_1(points, &indices, alloc);
    if colinear_len == indices.len() {
        return Convex3Result::Dim1(indices);
    }

    let coplanar_len = hull_2(points, &indices, colinear_len, alloc);

    let triangles = if coplanar_len == 3 {
        let mut triangles = Vec::with_capacity_in(3, alloc);
        triangles.extend(&indices[..3]);
        triangles
    } else {
        let axis = max_comp_in_tri_normal(
            point(points, indices[0]),
            point(points, indices[1]),
            point(points, indices[colinear_len]),
            alloc,
        );

        let mut points_2d = Vec::new_in(alloc);

        points_2d.extend(indices[..coplanar_len].into_iter().flat_map(|&idx| {
            let p = point(points, idx);
            [p[axis + 1], p[axis + 2]]
        }));
        let mut triangles = Vec::new_in(alloc);
        triangles.extend(
            triangulate(&points_2d, &[], alloc)
                .into_iter()
                .map(|idx| indices[idx]),
        );
        triangles
    };
    if colinear_len == indices.len() {
        return Convex3Result::Dim2(triangles);
    }
    Convex3Result::Dim3(hull_3(points, &indices, coplanar_len, triangles, alloc))
}

fn hull_1<A: Allocator + Copy>(points: &[f64], indices: &[usize], alloc: A) -> usize {
    let pa = point(points, indices[0]);
    let pb = point(points, indices[1]);
    for i in 2..indices.len() {
        let pc = point(points, indices[i]);
        if orient2d(pa, pb, pc, alloc) != 0.0 {
            return i;
        }
    }
    return indices.len();
}

fn hull_2<A: Allocator + Copy>(points: &[f64], indices: &[usize], start: usize, alloc: A) -> usize {
    let pa = point(points, indices[0]);
    let pb = point(points, indices[1]);
    let pc = point(points, indices[start]);
    for i in (start + 1)..indices.len() {
        let pd = point(points, indices[i]);
        if orient3d(pa, pb, pc, pd, alloc) != 0.0 {
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
) -> Vec<usize, A> {
    let is_pos = {
        let pa = point(points, triangles[0]);
        let pb = point(points, triangles[1]);
        let pc = point(points, triangles[2]);
        let pd = point(points, indices[start]);
        orient3d(pa, pb, pc, pd, alloc) > 0.0
    };

    if is_pos {
        triangles.reverse();
    }
    let mut mesh = ManifoldMesh::new(triangles.chunks(3));
    mesh.new_vertices(points.len() / 3 - mesh.n_vertices_capacity());

    let mut first_bot_hid = (mesh.n_halfedges_capacity() - 1).into();
    close_hull(&mut mesh, first_bot_hid, indices[start].into());

    let mut visited = Vec::new_in(alloc);
    visited.resize(mesh.n_faces_capacity(), false);
    let mut visible_faces = VecDeque::new_in(alloc);
    let mut removed_faces = Vec::new_in(alloc);
    let mut visited_faces = Vec::new_in(alloc);
    for &pid in &indices[(start + 1)..] {
        let pd = point(points, pid);
        let vid = pid.into();
        let visible = |fid: FaceId| {
            let he = mesh.face(fid).halfedge();

            let pa = point(points, he.to().0);
            let he = he.next();
            let pb = point(points, he.to().0);
            let he = he.next();
            let pc = point(points, he.to().0);
            orient3d(pa, pb, pc, pd, alloc) > 0.0
        };
        visible_faces.clear();
        removed_faces.clear();
        visited_faces.clear();
        visible_faces.extend(mesh.vertex(vid).incoming_halfedges().filter_map(|he| {
            let fid = mesh.he_face(*he);
            debug_assert!(fid.valid());
            visited[fid.0] = true;
            visited_faces.push(fid);
            if visible(fid) {
                removed_faces.push(fid);
                Some(fid)
            } else {
                None
            }
        }));

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
                visited_faces.push(twin_fid);

                if visible(twin_fid) {
                    visible_faces.push_back(twin_fid);
                    removed_faces.push(twin_fid);
                } else {
                    if !first_bot_hid.valid() {
                        first_bot_hid = *he;
                    }
                }
            }
        }

        debug_assert!(first_bot_hid.valid());
        for &fid in &removed_faces {
            mesh.remove_face(fid);
        }
        debug_assert!(mesh.he_is_boundary(first_bot_hid));

        for &fid in &visited_faces {
            visited[fid.0] = false;
        }

        close_hull(&mut mesh, first_bot_hid, vid);
    }

    Vec::new_in(alloc)
}

fn close_hull(mesh: &mut ManifoldMesh, first_hid: HalfedgeId, vid: VertexId) {
    debug_assert!(mesh.he_is_boundary(first_hid));

    let mut curr_bot_hid = first_hid;
    let mut first_side_hid = HalfedgeId::default();
    let mut prev_side_hid = HalfedgeId::default();
    loop {
        let next_bot_hid = mesh.he_next(curr_bot_hid);
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
