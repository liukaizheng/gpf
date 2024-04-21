use std::alloc::Allocator;

use hashbrown::HashMap;

use bumpalo::Bump;

use itertools::Itertools;

use crate::{
    predicates::{
        self, double_to_sign, inner_segment_cross_inner_triangle, inner_segment_cross_triangle,
        inner_segments_cross, max_comp_in_tri_normal, point_in_inner_segment, point_in_triangle,
        same_half_plane, Orientation,
    },
    triangle::{TetMesh, TriFace, VPIVOT},
    INVALID_IND,
};

pub(crate) struct Constraints {
    pub(crate) triangles: Vec<usize>,
    pub(crate) n_ori_triangles: usize,
}

impl Constraints {
    pub fn new(triangles: Vec<usize>) -> Self {
        let n_ori_triangles = triangles.len() / 3;
        Self {
            triangles,
            n_ori_triangles,
        }
    }

    #[inline(always)]
    pub(crate) fn triangle(&self, idx: usize) -> &[usize] {
        let start = idx * 3;
        &self.triangles[start..(start + 3)]
    }

    pub fn place_virtual_constraints(&mut self, tet_mesh: &TetMesh) {
        // make edge map: (va, vb) -> (tid, hid)
        let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for (idx, tri) in self.triangles.chunks(3).enumerate() {
            for (i, (&va, &vb)) in tri.iter().circular_tuple_windows().enumerate() {
                let key = if va < vb { (va, vb) } else { (vb, va) };
                if let Some(vec) = edge_map.get_mut(&key) {
                    vec.push((idx, i));
                } else {
                    edge_map.insert(key, vec![(idx, i)]);
                }
            }
        }
        let mut bump = Bump::new();
        let mut edges = edge_map.values().map(|v| v.clone()).collect::<Vec<_>>();
        for edge in &mut edges {
            edge.sort();
        }
        edges.sort();
        for halfedges in edges {
            bump.reset();
            self.add_virtual_constraint(tet_mesh, halfedges, &bump);
        }
    }

    fn add_virtual_constraint(
        &mut self,
        tet_mesh: &TetMesh,
        halfedges: Vec<(usize, usize)>,
        bump: &Bump,
    ) {
        // At least one half edge has been guaranteed
        let mut apex: usize;
        {
            let &(tid0, hid0) = &halfedges[0];
            let triangle0 = self.triangle(tid0);
            apex = triangle0[hid0];
            for &(tid, _) in &halfedges[1..] {
                let tri = self.triangle(tid);
                // unplannar: return
                if tet_mesh.orient3d(apex, tri[0], tri[1], tri[2], bump) != 0.0 {
                    return;
                }
            }

            // all triangles on the same plane
            if halfedges.len() > 1 {
                let axis = max_comp_in_tri_normal(
                    tet_mesh.point(triangle0[0]),
                    tet_mesh.point(triangle0[1]),
                    tet_mesh.point(triangle0[2]),
                    bump,
                );
                fn point2(point: &[f64], axis: usize) -> [f64; 2] {
                    [point[(axis + 1) % 3], point[(axis + 2) % 3]]
                }
                let pa = point2(tet_mesh.point(triangle0[(hid0 + 1) % 3]), axis);
                let pb = point2(tet_mesh.point(triangle0[(hid0 + 2) % 3]), axis);
                let mut pc = point2(tet_mesh.point(apex), axis);
                let base_ori = double_to_sign(predicates::orient2d(&pa, &pb, &pc, bump));
                for &(tid, hid) in &halfedges[1..] {
                    let triangle = self.triangle(tid);
                    apex = triangle[hid];
                    pc = point2(tet_mesh.point(apex), axis);
                    if double_to_sign(predicates::orient2d(&pa, &pb, &pc, bump)) != base_ori {
                        return;
                    }
                }
            }

            // all triangles are on the same side and coplanar
            let tet_id = tet_mesh.p2t[triangle0[0]];
            let tet = &tet_mesh.tets[tet_id].data;
            apex = tet[0];
            // if the first three times failed, then the next time must succeed
            for i in 1..4 {
                if apex == triangle0[0]
                    || apex == triangle0[1]
                    || apex == triangle0[2]
                    || tet_mesh.orient3d(apex, triangle0[0], triangle0[1], triangle0[2], bump)
                        == 0.0
                {
                    apex = tet[i];
                } else {
                    break;
                }
            }
        }
        let (tid, hid) = halfedges[0];
        let (org, dest) = {
            let tri = self.triangle(tid);
            (tri[(hid + 1) % 3], tri[(hid + 2) % 3])
        };
        self.triangles.push(org);
        self.triangles.push(dest);
        self.triangles.push(apex);
    }

    pub fn insert_constraints<'a>(&self, mesh: &mut TetMesh) -> [Vec<Vec<usize>>; 5] {
        let mut tet_marks = [
            vec![Vec::<usize>::new(); mesh.tets.len()],
            vec![Vec::<usize>::new(); mesh.tets.len()],
            vec![Vec::<usize>::new(); mesh.tets.len()],
            vec![Vec::<usize>::new(); mesh.tets.len()],
            vec![Vec::<usize>::new(); mesh.tets.len()],
        ];
        let mut bump = Bump::new();
        for (i, triangle) in self.triangles.chunks(3).enumerate() {
            bump.reset();
            if i == 1297 {
                println!("1297");
            }
            let tet_face = triangle_at_tet(mesh, triangle, &bump);
            if tet_face.tet != INVALID_IND {
                tet_marks[tet_face.ver & 3][tet_face.tet].push(i);
                let nei = &mesh.tets[tet_face.tet].nei[tet_face.ver & 3];
                tet_marks[nei.ver & 3][nei.tet].push(i);
                continue;
            }
            let mut intersect_info = IntersectInfo {
                intersected: Vec::new_in(&bump),
                visited: std::vec::from_elem_in(false, mesh.tets.len(), &bump),
            };
            constraint_sides_intersections(mesh, triangle, &mut intersect_info);
            set_improper_intersections(mesh, i, triangle, &mut tet_marks, &mut intersect_info);
            interior_intersections(mesh, i, triangle, &mut tet_marks, &mut intersect_info);
        }
        tet_marks
    }
}

struct IntersectInfo<A: Allocator + Copy> {
    intersected: Vec<usize, A>,
    visited: Vec<bool, A>,
}

#[inline(always)]
fn vert_inner_segment_cross_inner_triangle<A: Allocator + Copy>(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    v1: usize,
    v2: usize,
    v3: usize,
    allocator: A,
) -> bool {
    if u1 == v1 || u1 == v2 || u1 == v3 || u2 == v1 || u2 == v2 || u2 == v3 {
        return false;
    }
    return inner_segment_cross_inner_triangle(
        mesh.point(u1),
        mesh.point(u2),
        mesh.point(v1),
        mesh.point(v2),
        mesh.point(v3),
        allocator,
    );
}

#[inline(always)]
fn vert_inner_segment_cross_triangle<A: Allocator + Copy>(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    tri: &[usize],
    allocator: A,
) -> bool {
    if u1 == tri[0] || u1 == tri[1] || u1 == tri[2] || u2 == tri[0] || u2 == tri[1] || u2 == tri[2]
    {
        return false;
    }
    return inner_segment_cross_triangle(
        mesh.point(u1),
        mesh.point(u2),
        mesh.point(tri[0]),
        mesh.point(tri[1]),
        mesh.point(tri[2]),
        allocator,
    );
}

#[inline(always)]
fn vert_inner_segments_cross<A: Allocator + Copy>(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    v1: usize,
    v2: usize,
    allocator: A,
) -> bool {
    if u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2 {
        return false;
    }
    return inner_segments_cross(
        mesh.point(u1),
        mesh.point(u2),
        mesh.point(v1),
        mesh.point(v2),
        allocator,
    );
}

#[inline(always)]
fn verts_in_same_half_space<A: Allocator + Copy>(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    v1: usize,
    v2: usize,
    v3: usize,
    allocator: A,
) -> bool {
    double_to_sign(mesh.orient3d(u1, v1, v2, v3, allocator))
        == double_to_sign(mesh.orient3d(u2, v1, v2, v3, allocator))
}

#[inline(always)]
fn verts_same_half_plane<A: Allocator + Copy>(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    v1: usize,
    v2: usize,
    allocator: A,
) -> bool {
    if u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2 {
        return false;
    }
    return same_half_plane(
        mesh.point(u1),
        mesh.point(u2),
        mesh.point(v1),
        mesh.point(v2),
        allocator,
    );
}

#[inline(always)]
fn vert_point_in_inner_segment<A: Allocator + Copy>(
    mesh: &TetMesh,
    u: usize,
    v1: usize,
    v2: usize,
    allocator: A,
) -> bool {
    return u != v1
        && u != v2
        && point_in_inner_segment(mesh.point(u), mesh.point(v1), mesh.point(v2), allocator);
}
#[inline(always)]
fn vert_point_in_segment<A: Allocator + Copy>(
    mesh: &TetMesh,
    u: usize,
    v1: usize,
    v2: usize,
    allocator: A,
) -> bool {
    return u == v1
        || u == v2
        || point_in_inner_segment(mesh.point(u), mesh.point(v1), mesh.point(v2), allocator);
}

#[inline(always)]
fn vert_point_in_inner_triangle<A: Allocator + Copy>(
    mesh: &TetMesh,
    p: usize,
    tri: &[usize],
    allocator: A,
) -> bool {
    return p != tri[0]
        && p != tri[1]
        && p != tri[2]
        && point_in_triangle(
            mesh.point(p),
            mesh.point(tri[0]),
            mesh.point(tri[1]),
            mesh.point(tri[2]),
            allocator,
        );
}

#[inline(always)]
fn vert_point_in_triangle<A: Allocator + Copy>(
    mesh: &TetMesh,
    p: usize,
    tri: &[usize],
    allocator: A,
) -> bool {
    return p == tri[0]
        || p == tri[1]
        || p == tri[2]
        || point_in_triangle(
            mesh.point(p),
            mesh.point(tri[0]),
            mesh.point(tri[1]),
            mesh.point(tri[2]),
            allocator,
        );
}

fn triangle_at_tet(mesh: &mut TetMesh, tri: &[usize], bump: &Bump) -> TriFace {
    let mut tets = [mesh.p2t[tri[0]]].to_vec_in(bump);
    mesh.mark_test(tets[0]);

    // find the second point in triangle
    let mut pos3 = 3; // position of the third point in tet
    let mut edge = TriFace::default();
    let mut idx = 0;
    while idx < tets.len() {
        let tid = tets[idx];
        let pos1 = mesh.tets[tid].index(tri[0]).unwrap();
        let mut cur = TriFace::new(tid, VPIVOT[pos1]);
        for _ in 0..3 {
            let dest = mesh.dest(&cur);
            if dest == tri[1] || dest == tri[2] {
                pos3 = if dest == tri[1] { 2 } else { 1 };
                edge.set(cur.tet, cur.ver);
                break;
            } else {
                cur.eprev_esym_self();
            }
        }
        if pos3 != 3 {
            break;
        }

        for i in 0..4 {
            if i == pos1 {
                break;
            }
            let nei = mesh.tets[tid].nei[i].tet;
            if !mesh.is_hull_tet(nei) && !mesh.mark_tested(nei) {
                mesh.mark_test(nei);
                tets.push(nei);
            }
        }
        idx += 1;
    }

    for tid in tets {
        mesh.unmark_test(tid);
    }

    if pos3 == 3 {
        return TriFace::default();
    }

    // find the third point in triangle
    let mut spin = edge.clone();
    loop {
        if mesh.apex(&spin) == tri[pos3] {
            return spin;
        }
        mesh.fnext_self(&mut spin);
        if spin.tet == edge.tet {
            break;
        }
    }
    TriFace::default()
}

/// find all tets which intersect with triangle side
fn constraint_sides_intersections<A: Allocator + Copy>(
    mesh: &mut TetMesh,
    triangle: &[usize],
    info: &mut IntersectInfo<A>,
) {
    for (&va, &vb) in triangle.iter().circular_tuple_windows() {
        let (mut tid, mut connect_verts) = vert_on_constraint_side(mesh, va, vb, INVALID_IND, info);
        while connect_verts[0] != vb {
            match connect_verts.len() {
                1 => {
                    (tid, connect_verts) =
                        vert_on_constraint_side(mesh, connect_verts[0], vb, tid, info);
                }
                2 => {
                    let mut edge =
                        TriFace::new(tid, VPIVOT[mesh.tets[tid].index(connect_verts[0]).unwrap()]);
                    // the third must succeed
                    for _ in 0..2 {
                        if mesh.dest(&edge) == connect_verts[1] {
                            break;
                        }
                        edge.eprev_esym_self();
                    }
                    (tid, connect_verts) =
                        edge_cross_constraint_side(mesh, edge, va, vb, tid, info);
                }
                3 => {
                    let mut edge =
                        TriFace::new(tid, VPIVOT[mesh.tets[tid].index(connect_verts[0]).unwrap()]);
                    for _ in 0..2 {
                        let (dest, apex) = (mesh.dest(&edge), mesh.apex(&edge));
                        if (dest == connect_verts[1] && apex == connect_verts[2])
                            || (dest == connect_verts[2] && apex == connect_verts[1])
                        {
                            break;
                        }
                        edge.eprev_esym_self();
                    }
                    (tid, connect_verts) =
                        triangle_pierced_by_constraint_side(mesh, edge, va, vb, info);
                }
                _ => unreachable!("never"),
            }
        }
    }
}

fn vert_on_constraint_side<A: Allocator + Copy>(
    mesh: &mut TetMesh,
    va: usize,
    vb: usize,
    prev_tid: usize,
    info: &mut IntersectInfo<A>,
) -> (usize, Vec<usize, A>) {
    let bump = *info.intersected.allocator();
    let inc_tets = mesh.incident(va, bump);
    for &tid in &inc_tets {
        if !info.visited[tid] {
            info.intersected.push(tid);
            info.visited[tid] = true;
        }
    }

    for tid in inc_tets {
        if tid == prev_tid {
            continue;
        }
        let tet = &mesh.tets[tid];
        if tet.index(vb).is_some() {
            return (tid, [vb].to_vec_in(bump));
        }
        // the vertices of opposite face aganist "va"
        let mut oppo_verts = Vec::with_capacity_in(3, bump);
        oppo_verts.extend(
            tet.data
                .iter()
                .filter_map(|&vid| if vid != va { Some(vid) } else { None }),
        );

        if vert_inner_segment_cross_inner_triangle(
            mesh,
            va,
            vb,
            oppo_verts[0],
            oppo_verts[1],
            oppo_verts[2],
            bump,
        ) {
            return (tid, oppo_verts);
        }

        for (&ua, &ub) in oppo_verts.iter().circular_tuple_windows() {
            if vert_inner_segments_cross(mesh, va, vb, ua, ub, bump) {
                return (tid, [ua, ub].to_vec_in(bump));
            }
        }

        for vid in oppo_verts {
            if vert_point_in_inner_segment(mesh, vid, va, vb, bump) {
                return (tid, [vid].to_vec_in(bump));
            }
        }
    }
    unreachable!("cannot find intersect tet by line segment");
}

fn edge_cross_constraint_side<A: Allocator + Copy>(
    mesh: &mut TetMesh,
    edge: TriFace,
    ea: usize,
    eb: usize,
    prev_tid: usize,
    info: &mut IntersectInfo<A>,
) -> (usize, Vec<usize, A>) {
    let bump = *info.visited.allocator();
    let mut spin = edge.clone();
    mesh.fnext_self(&mut spin);
    while spin.tet != edge.tet {
        let tid = spin.tet;
        if !mesh.is_hull_tet(tid) {
            if !info.visited[tid] {
                info.visited[tid] = true;
                info.intersected.push(tid);
            }
        }
        mesh.fnext_self(&mut spin);
    }

    // now spin == edge
    mesh.fnext_self(&mut spin);
    let (va, vb) = (mesh.org(&edge), mesh.dest(&edge));
    while spin.tet != edge.tet {
        let tid = spin.tet;

        if tid != prev_tid && !mesh.is_hull_tet(tid) {
            let (vc, vd) = (mesh.apex(&spin), mesh.oppo(&spin));

            // segment is one of edges of tet
            if vc == eb || vd == eb {
                return (tid, [eb].to_vec_in(bump));
            }

            // segment and triangle properly intersect
            if vert_inner_segment_cross_inner_triangle(mesh, ea, eb, vc, vd, va, bump)
                && verts_in_same_half_space(mesh, vb, ea, vc, vd, va, bump)
            {
                return (tid, [vc, vd, va].to_vec_in(bump));
            }
            if vert_inner_segment_cross_inner_triangle(mesh, ea, eb, vc, vd, vb, bump)
                && verts_in_same_half_space(mesh, va, ea, vc, vd, vb, bump)
            {
                return (tid, [vc, vd, vb].to_vec_in(bump));
            }

            // segment and segment properly intersect
            if vert_inner_segments_cross(mesh, ea, eb, vc, vd, bump)
                && verts_in_same_half_space(mesh, ea, va, vc, vd, vb, bump)
            {
                return (tid, [vc, vd].to_vec_in(bump));
            }

            // tet vertex on segment
            if vert_point_in_inner_segment(mesh, vc, ea, eb, bump)
                && verts_same_half_plane(mesh, vc, eb, va, vb, bump)
            {
                return (tid, [vc].to_vec_in(bump));
            }
            if vert_point_in_inner_segment(mesh, vd, ea, eb, bump)
                && verts_same_half_plane(mesh, vd, eb, va, vb, bump)
            {
                return (tid, [vd].to_vec_in(bump));
            }

            // segment and segment properly intersect
            for v1 in [va, vb] {
                for v2 in [vc, vd] {
                    if vert_inner_segments_cross(mesh, ea, eb, v1, v2, bump)
                        && verts_same_half_plane(mesh, v2, eb, va, vb, bump)
                    {
                        return (tid, [v1, v2].to_vec_in(bump));
                    }
                }
            }
        }
        mesh.fnext_self(&mut spin);
    }

    panic!("We cannot find the next tet when cross some tet edge")
}
fn triangle_pierced_by_constraint_side<A: Allocator + Copy>(
    mesh: &mut TetMesh,
    mut tri: TriFace,
    ea: usize,
    eb: usize,
    info: &mut IntersectInfo<A>,
) -> (usize, Vec<usize, A>) {
    let next_tet = &mesh.tets[tri.tet].nei[tri.ver & 3];
    let next_tid = next_tet.tet;
    if !info.visited[next_tid] {
        info.visited[next_tid] = true;
        info.intersected.push(next_tid);
    }

    let v_oppo = mesh.oppo(next_tet);
    let bump = *info.visited.allocator();
    if v_oppo == eb {
        return (next_tid, [eb].to_vec_in(bump));
    }

    for _ in 0..3 {
        let va = mesh.org(&tri);
        let vb = mesh.dest(&tri);
        if vert_inner_segment_cross_inner_triangle(mesh, ea, eb, va, vb, v_oppo, bump) {
            return (next_tid, [va, vb, v_oppo].to_vec_in(bump));
        }
        tri.enext_self();
    }

    for v in [mesh.org(&tri), mesh.dest(&tri), mesh.apex(&tri)] {
        if vert_inner_segments_cross(mesh, ea, eb, v, v_oppo, bump) {
            return (next_tid, [v, v_oppo].to_vec_in(bump));
        }
    }

    // it must be that [ea, eb] pass through v_oppo
    return (next_tid, [v_oppo].to_vec_in(bump));
}

fn set_improper_intersections<A: Allocator + Copy>(
    mesh: &TetMesh,
    cid: usize,
    triangle: &[usize],
    tet_marks: &mut [Vec<Vec<usize>>],
    info: &mut IntersectInfo<A>,
) {
    let bump = *info.intersected.allocator();
    let lookup_ori = |vid: usize, vert_oris: &mut [f64], used_verts: &mut Vec<usize, A>| -> f64 {
        let mut ori = vert_oris[vid];
        if !ori.is_nan() {
            return ori;
        } else {
            ori = mesh.orient3d(triangle[0], triangle[1], triangle[2], vid, bump);
            vert_oris[vid] = ori;
            used_verts.push(vid);
            return ori;
        }
    };
    let mut vert_oris = std::vec::from_elem_in(f64::NAN, mesh.tets.len(), bump);
    let mut used_verts = Vec::new_in(bump);

    for &tid in &info.intersected {
        let tet = &mesh.tets[tid];

        let mut zero = Vec::new_in(bump);
        let mut pos = Vec::new_in(bump);
        let mut neg = Vec::new_in(bump);
        for &vid in &tet.data {
            let ori = lookup_ori(vid, &mut vert_oris, &mut used_verts);
            if ori == 0.0 {
                zero.push(vid);
            } else if ori > 0.0 {
                pos.push(vid);
            } else {
                neg.push(vid);
            }
        }

        match zero.len() {
            3 => {
                let oppo_index = tet
                    .index(if pos.is_empty() { neg[0] } else { pos[0] })
                    .unwrap();
                if constraint_intersect_face(mesh, triangle, &zero, bump) {
                    tet_marks[oppo_index][tid].push(cid);
                }
            }
            2 => {
                if neg.len() != 1 {
                    continue;
                }
                let non_planar_v = [neg[0], pos[0]];
                if constraint_intersect_edge(mesh, triangle, &zero, &non_planar_v, bump) {
                    tet_marks[4][tid].push(cid);
                }
            }
            1 => {
                if neg.is_empty() || pos.is_empty() {
                    continue;
                }
                let (under_vs, over_v) = if neg.len() > pos.len() {
                    (&neg, pos[0])
                } else {
                    (&pos, neg[0])
                };
                if constraint_intersect_vertex(mesh, triangle, zero[0], over_v, under_vs, bump) {
                    tet_marks[4][tid].push(cid);
                }
            }
            _ => {
                tet_marks[4][tid].push(cid);
            }
        }
    }
}

fn constraint_intersect_face<A: Allocator + Copy>(
    mesh: &TetMesh,
    c: &[usize],
    t: &[usize],
    bump: A,
) -> bool {
    if vert_point_in_triangle(mesh, t[0], c, bump)
        && vert_point_in_triangle(mesh, t[1], c, bump)
        && vert_point_in_triangle(mesh, t[2], c, bump)
    {
        return true;
    }

    if vert_inner_segments_cross(mesh, t[0], t[1], c[0], c[1], bump)
        || vert_inner_segments_cross(mesh, t[0], t[1], c[1], c[2], bump)
        || vert_inner_segments_cross(mesh, t[0], t[1], c[2], c[0], bump)
        || vert_inner_segments_cross(mesh, t[1], t[2], c[0], c[1], bump)
        || vert_inner_segments_cross(mesh, t[1], t[2], c[1], c[2], bump)
        || vert_inner_segments_cross(mesh, t[1], t[2], c[2], c[0], bump)
        || vert_inner_segments_cross(mesh, t[2], t[0], c[0], c[1], bump)
        || vert_inner_segments_cross(mesh, t[2], t[0], c[1], c[2], bump)
        || vert_inner_segments_cross(mesh, t[2], t[0], c[2], c[0], bump)
    {
        return true;
    };
    return false;
}

fn constraint_intersect_edge<A: Allocator + Copy>(
    mesh: &TetMesh,
    c: &[usize],
    coplanar_v: &[usize],
    non_planar_v: &[usize],
    bump: A,
) -> bool {
    let (mut va_in, mut vb_in) = (
        vert_point_in_inner_triangle(mesh, coplanar_v[0], c, bump),
        vert_point_in_inner_triangle(mesh, coplanar_v[1], c, bump),
    );
    // va and vb both in the interior of triangle, intersection must be improper.
    if va_in && vb_in {
        return true;
    }
    let va_on_edge = vert_point_in_segment(mesh, coplanar_v[0], c[1], c[2], bump) as u8
        | ((vert_point_in_segment(mesh, coplanar_v[0], c[2], c[0], bump) as u8) << 1)
        | ((vert_point_in_segment(mesh, coplanar_v[0], c[0], c[1], bump) as u8) << 2);
    let vb_on_edge = vert_point_in_segment(mesh, coplanar_v[1], c[1], c[2], bump) as u8
        | ((vert_point_in_segment(mesh, coplanar_v[1], c[2], c[0], bump) as u8) << 1)
        | ((vert_point_in_segment(mesh, coplanar_v[1], c[0], c[1], bump) as u8) << 2);
    va_in |= va_on_edge > 0;
    vb_in |= vb_on_edge > 0;
    if va_in && vb_in {
        let on_same_edge = va_on_edge & vb_on_edge;
        // va and vb are not on same edge
        if on_same_edge == 0 {
            return true;
        }
        // edge [vc vd] pierces triangle
        if vert_inner_segment_cross_triangle(mesh, non_planar_v[0], non_planar_v[1], c, bump) {
            return true;
        }
        let bit = on_same_edge.trailing_zeros() as usize;
        let ea = c[(bit + 1) % 3];
        let eb = c[(bit + 2) % 3];
        for ev in [ea, eb] {
            for &tv in coplanar_v {
                if vert_inner_segment_cross_inner_triangle(
                    mesh,
                    ev,
                    c[bit],
                    non_planar_v[0],
                    non_planar_v[1],
                    tv,
                    bump,
                ) {
                    return true;
                }
            }
        }
        // otherwise.. the intersection is proper
        return false;
    } else if va_in || vb_in {
        // edge [va vb] crosses triangle edge
        if vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[0], c[1], bump)
            || vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[1], c[2], bump)
            || vert_inner_segments_cross(mesh, coplanar_v[0], coplanar_v[1], c[2], c[0], bump)
        {
            return true;
        }
        // edge [vc vd] pierces triangle
        if vert_inner_segment_cross_triangle(mesh, non_planar_v[0], non_planar_v[1], c, bump) {
            return true;
        }
        for (&ea, &eb) in c.iter().circular_tuple_windows() {
            for &tv in coplanar_v {
                if vert_inner_segment_cross_inner_triangle(
                    mesh,
                    ea,
                    eb,
                    non_planar_v[0],
                    non_planar_v[1],
                    tv,
                    bump,
                ) {
                    return true;
                }
            }
        }
        // otherwise.. the intersection is proper
        return false;
    } else {
        // otherwise.. since an intersection exists, it must be improper
        return true;
    }
}

fn constraint_intersect_vertex<A: Allocator + Copy>(
    mesh: &TetMesh,
    c: &[usize],
    tv: usize,
    over_v: usize,
    under_vs: &[usize],
    bump: A,
) -> bool {
    if vert_point_in_inner_triangle(mesh, tv, c, bump) {
        return true;
    }
    let vert_on_segments = vert_point_in_segment(mesh, tv, c[0], c[1], bump)
        || vert_point_in_segment(mesh, tv, c[1], c[2], bump)
        || vert_point_in_segment(mesh, tv, c[2], c[0], bump);
    // tv isn't on triangle, but intersection exist, so it's improper
    if !vert_on_segments {
        return true;
    }
    // now tv is on a boundary of triangle
    if vert_inner_segment_cross_triangle(mesh, over_v, under_vs[0], c, bump)
        && vert_inner_segment_cross_triangle(mesh, over_v, under_vs[1], c, bump)
    {
        return true;
    }
    for (&ca, &cb) in c.iter().circular_tuple_windows() {
        if vert_inner_segment_cross_inner_triangle(mesh, ca, cb, tv, over_v, under_vs[0], bump)
            || vert_inner_segment_cross_inner_triangle(mesh, ca, cb, tv, over_v, under_vs[1], bump)
            || vert_inner_segment_cross_inner_triangle(
                mesh,
                ca,
                cb,
                over_v,
                under_vs[0],
                under_vs[1],
                bump,
            )
        {
            return true;
        }
    }
    // otherwise... proper intersection
    return false;
}

fn interior_intersections<A: Allocator + Copy>(
    mesh: &TetMesh,
    cid: usize,
    triangle: &[usize],
    tet_marks: &mut [Vec<Vec<usize>>],
    info: &mut IntersectInfo<A>,
) {
    let mut idx = 0;
    let bump = *info.visited.allocator();
    while idx < info.intersected.len() {
        let tid = info.intersected[idx];
        for adj_tid in mesh.tets[tid].nei.iter().map(|f| f.tet) {
            if !info.visited[adj_tid] {
                info.visited[adj_tid] = true;
                let t = constraint_intersect_tet_type(mesh, triangle, adj_tid, bump);
                if t < 5 {
                    tet_marks[t][adj_tid].push(cid);
                    info.intersected.push(adj_tid);
                }
            }
        }
        idx += 1;
    }
}

fn constraint_intersect_tet_type<A: Allocator + Copy>(
    mesh: &TetMesh,
    c: &[usize],
    tid: usize,
    bump: A,
) -> usize {
    if mesh.is_hull_tet(tid) {
        return 5;
    }
    let verts = &mesh.tets[tid].data;
    let mut oris = [Orientation::Undefined; 4];
    let mut in_tris = [false; 4];
    for i in 0..4 {
        oris[i] = double_to_sign(mesh.orient3d(verts[i], c[0], c[1], c[2], bump));
        in_tris[i] = oris[i] == Orientation::Zero && vert_point_in_inner_triangle(mesh, verts[i], c, bump);
    }

    let mut indices = [0, 1, 2, 3];
    let mid = indices.iter_mut().partition_in_place(|&i| in_tris[i]); // indices[..mid]: in triangle.
    let (verts_in_c, verts_out_c) = indices.split_at(mid);
    match verts_in_c.len() {
        3 => {
            return verts_out_c[0];
        }
        2 => {
            if oris[verts_out_c[0]] == oris[verts_out_c[1]] {
                return 5;
            } else {
                return 4;
            }
        }
        1 => {
            let ori1 = oris[verts_out_c[0]];
            let ori2 = oris[verts_out_c[1]];
            let ori3 = oris[verts_out_c[2]];
            if ori1 == ori2 && ori1 == ori3 {
                return 5;
            } else {
                return 4;
            }
        }
        _ => {
            for i in 0..3 {
                for j in (i + 1)..4 {
                    if oris[i] != oris[j]
                        && vert_inner_segment_cross_inner_triangle(
                            mesh, verts[i], verts[j], c[0], c[1], c[2], bump,
                        )
                    {
                        return 4;
                    }
                }
            }
            return 5;
        }
    }
}
