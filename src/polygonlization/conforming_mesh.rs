use std::collections::HashMap;

use bumpalo::collections::Vec;

use crate::{
    predicates::{self, double_to_sign, max_comp_in_tri_normal},
    triangle::TetMesh,
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

    fn triangle(&self, idx: usize) -> &[usize] {
        let start = idx * 3;
        &self.triangles[start..(start + 3)]
    }

    pub fn place_virtual_constraints<'a, 'bb: 'a>(&mut self, tet_mesh: &TetMesh<'a, 'bb>) {
        // make edge map: (va, vb) -> (tid, hid)
        let bump = self.triangles.bump();
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
            self.add_virtual_constrait(tet_mesh, halfedges);
        }
    }

    fn add_virtual_constrait<'a, 'bb: 'a>(
        &mut self,
        tet_mesh: &TetMesh<'a, 'bb>,
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
                let bump = self.triangles.bump();
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

            // all triangles are on the same side and coplannar
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
}
