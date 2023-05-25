use std::collections::HashMap;

use bumpalo::{collections::Vec, Bump};

use crate::{
    predicates::{
        self, double_to_sign, inner_segment_cross_inner_triangle, inner_segment_cross_triangle,
        max_comp_in_tri_normal,
    },
    triangle::{TetMesh, TriFace, VPIVOT},
    INVALID_IND,
};
pub struct Constraints<'b> {
    triangles: Vec<'b, usize>,
    n_ori_triangles: usize,
}

impl<'b> Constraints<'b> {
    pub fn new(triangles: Vec<'b, usize>) -> Self {
        let n_ori_triangles = triangles.len() / 3;
        Self {
            triangles,
            n_ori_triangles,
        }
    }

    #[inline(always)]
    fn triangle(&self, idx: usize) -> &[usize] {
        let start = idx * 3;
        &self.triangles[start..(start + 3)]
    }

    #[inline(always)]
    fn bump(&self) -> &'b Bump {
        self.triangles.bump()
    }

    pub fn place_virtual_constraints<'a>(&mut self, tet_mesh: &TetMesh<'a, 'b>) {
        // make edge map: (va, vb) -> (tid, hid)
        let bump = self.bump();
        let mut edge_map: HashMap<(usize, usize), Vec<'b, (usize, usize)>> = HashMap::new();
        for (idx, tri) in self.triangles.chunks(3).enumerate() {
            for i in 0..3 {
                let va = tri[(i + 1) % 3];
                let vb = tri[(i + 2) % 3];
                let key = if va < vb { (va, vb) } else { (vb, va) };
                if let Some(vec) = edge_map.get_mut(&key) {
                    vec.push((idx, i));
                } else {
                    edge_map.insert(key, bumpalo::vec![in bump;(idx, i)]);
                }
            }
        }
        for (_, halfedges) in edge_map {
            self.add_virtual_constraint(tet_mesh, halfedges);
        }
    }

    fn add_virtual_constraint<'a>(
        &mut self,
        tet_mesh: &TetMesh<'a, 'b>,
        halfedges: Vec<'b, (usize, usize)>,
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
                if tet_mesh.orient3d(apex, tri[0], tri[1], tri[2]) != 0.0 {
                    return;
                }
            }

            // all triangles on the same plane
            if halfedges.len() > 1 {
                let bump = self.bump();
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
                    || tet_mesh.orient3d(apex, triangle0[0], triangle0[1], triangle0[2]) == 0.0
                {
                    apex = tet[i + 1];
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

    pub fn insert_constraints<'a>(&self, mesh: &mut TetMesh<'a, 'b>) -> [Vec<Vec<usize>>; 5] {
        let bump = mesh.tets.bump();
        let mut tet_map = [
            bumpalo::vec![in bump; Vec::<'b, usize>::new_in(bump); mesh.tets.len()],
            bumpalo::vec![in bump; Vec::<'b, usize>::new_in(bump); mesh.tets.len()],
            bumpalo::vec![in bump; Vec::<'b, usize>::new_in(bump); mesh.tets.len()],
            bumpalo::vec![in bump; Vec::<'b, usize>::new_in(bump); mesh.tets.len()],
            bumpalo::vec![in bump; Vec::<'b, usize>::new_in(bump); mesh.tets.len()],
        ];
        for (i, triangle) in self.triangles.chunks(3).enumerate() {
            let tet_face = triangle_at_tet(mesh, triangle);
            if tet_face.tet != INVALID_IND {
                tet_map[tet_face.ver & 3][tet_face.tet].push(i);
                let nei = &mesh.tets[tet_face.tet].nei[tet_face.ver & 3];
                tet_map[nei.ver & 3][nei.tet].push(i);
            }
            let intersect_info = IntersectInfo {
                intersected: Vec::new_in(bump),
                visited: bumpalo::vec![in bump; false; mesh.tets.len()],
            };
        }
        tet_map
    }
}

struct IntersectInfo<'b> {
    intersected: Vec<'b, usize>,
    visited: Vec<'b, bool>,
}

#[inline(always)]
fn vert_inner_segment_cross_inner_triangle(
    mesh: &TetMesh,
    u1: usize,
    u2: usize,
    v1: usize,
    v2: usize,
    v3: usize,
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
        mesh.tets.bump(),
    );
}

#[inline(always)]
fn vert_inner_segment_cross_triangle(mesh: &TetMesh, u1: usize, u2: usize, tri: &[usize]) -> bool {
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
        mesh.tets.bump(),
    );
}

fn triangle_at_tet<'a, 'b: 'a>(mesh: &mut TetMesh<'a, 'b>, tri: &[usize]) -> TriFace {
    let bump = mesh.tets.bump();
    let mut tets = bumpalo::vec![in bump; mesh.p2t[0]];
    mesh.mark_test(tets[0]);

    // find the second point in triangle
    let mut pos3 = 3; // position of the third point in tet
    let mut edge = TriFace::default();
    let idx = 0;
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
            if mesh.is_hull_tet(nei) && !mesh.mark_tested(nei) {
                tets.push(nei);
            }
        }
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
fn constraint_sides_intersection(mesh: &mut TetMesh, triangle: &[usize]) {}

fn constraint_side_intersection_start<'b>(
    mesh: &mut TetMesh<'_, 'b>,
    va: usize,
    vb: usize,
    info: &mut IntersectInfo<'b>,
) -> (usize, Vec<'b, usize>) {
    let bump = mesh.tets.bump();
    let inc_tets = mesh.incident(va);
    for &tid in &inc_tets {
        if !info.visited[tid] {
            info.intersected.push(tid);
            info.visited[tid] = true;
        }
    }

    for tid in inc_tets {
        let tet = &mesh.tets[tid];
        if tet.index(vb).is_some() {
            return (tid, bumpalo::vec![in bump; vb]);
        }
        // the vertices of opposite face aganist "va"
        let oppo_verts = Vec::from_iter_in(
            tet.data
                .iter()
                .filter_map(|&vid| if vid != vb { Some(vid) } else { None }),
            bump,
        );

        if vert_inner_segment_cross_inner_triangle(
            mesh,
            va,
            vb,
            oppo_verts[0],
            oppo_verts[1],
            oppo_verts[2],
        ) {
            return (tid, oppo_verts);
        }
    }
    unreachable!("cannot find intersect tet by line segment");
}
