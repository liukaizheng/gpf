use core::panic;
use std::collections::{BinaryHeap, HashMap};

use bumpalo::Bump;
use itertools::Itertools;

use crate::{
    geometry::{Surf, Surface},
    math::{cross, cross_in, dot, square_norm, sub_short},
    mesh::{EdgeId, Mesh, SurfaceMesh},
    point,
    predicates::{double_to_sign, orient2d, orient3d, Orientation},
};

type BVec<'a, T> = bumpalo::collections::Vec<'a, T>;

use super::TetSet;

#[derive(Debug, Clone)]
struct EdgeAndLen {
    eid: EdgeId,
    len: f64,
}

impl PartialEq for EdgeAndLen {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.eid == other.eid && self.len == other.len
    }
}

impl Eq for EdgeAndLen {}

impl PartialOrd for EdgeAndLen {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.len.partial_cmp(&other.len)
    }
}

impl Ord for EdgeAndLen {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

struct SubdivisionData<'a> {
    surfaces: Vec<&'a Surf>,
    vals_and_grads: Vec<HashMap<usize, [f64; 4]>>,
    active_surfs: Vec<Vec<usize>>,
    queue: BinaryHeap<EdgeAndLen>,
}

pub(super) fn adaptive_subdivide(tets: &mut TetSet, surfaces: Vec<&Surf>, sq_eps: f64) {
    let mut check_bump = Bump::new();
    let mut vals_and_grads = vec![HashMap::<usize, [f64; 4]>::new(); surfaces.len()];
    // TODO: parallelize by rayon
    for (i, &surf) in surfaces.iter().enumerate() {
        for (vid, p) in tets.points.chunks(3).enumerate() {
            vals_and_grads[i].insert(vid, surf.eval(p));
        }
    }
    let active_surfs = (0..tets.tet_faces.len())
        .map(|_| Vec::from_iter(0..surfaces.len()))
        .collect_vec();
    let mut data = SubdivisionData {
        surfaces,
        vals_and_grads,
        active_surfs,
        queue: BinaryHeap::new(),
    };

    for tid in 0..tets.tet_faces.len() {
        check_bump.reset();
        push_longest_edge(tid, tets, &mut data, sq_eps, &check_bump);
    }

    let mut split_bump = Bump::new();
    while !data.queue.is_empty() {
        let EdgeAndLen { eid, len } = data.queue.pop().unwrap();
        if tets.square_edge_lengths[eid] != len {
            // edge changed
            continue;
        }
        split_bump.reset();
        let (new_vert, tet_pairs) = tets.split_edge(eid, &split_bump);
        for &[old_tet, _] in &tet_pairs {
            for &sid in &data.active_surfs[old_tet] {
                let vals_grads = &mut data.vals_and_grads[sid];
                if vals_grads.contains_key(&new_vert) {
                    continue;
                }
                vals_grads.insert(
                    new_vert.0,
                    data.surfaces[sid].eval(point(&tets.points, new_vert.0)),
                );
            }
            data.active_surfs.push(data.active_surfs[old_tet].clone());
        }
        for tid in tet_pairs.into_iter().flatten() {
            check_bump.reset();
            push_longest_edge(tid, tets, &mut data, sq_eps, &check_bump);
        }
    }
    println!("n tets: {}", tets.tet_faces.len());
}

fn push_longest_edge(
    tid: usize,
    tets: &mut TetSet,
    data: &mut SubdivisionData,
    sq_eps: f64,
    bump: &Bump,
) {
    if subdividable(tid, tets, data, sq_eps, bump) {
        let longest_eid = *tets.tet_edges[tid]
            .iter()
            .max_by(|&&ea, &&eb| {
                tets.square_edge_lengths[ea]
                    .partial_cmp(&tets.square_edge_lengths[eb])
                    .unwrap()
            })
            .unwrap();
        data.queue.push(EdgeAndLen {
            eid: longest_eid,
            len: tets.square_edge_lengths[longest_eid],
        });
    }
}

const C: [[f64; 4]; 16] = [
    [2.0 / 3.0, 1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 1.0 / 3.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0, 1.0 / 3.0],
    [0.0, 2.0 / 3.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0],
    [0.0, 0.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 2.0 / 3.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 0.0, 2.0 / 3.0],
    [0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
];

fn subdividable(
    tid: usize,
    tets: &mut TetSet,
    data: &mut SubdivisionData,
    sq_eps: f64,
    bump: &Bump,
) -> bool {
    let surfs = &mut data.active_surfs[tid];
    let mut active = Vec::with_capacity_in(surfs.len(), bump);
    active.resize(surfs.len(), true);
    let verts = tets.vertices_in(tid, bump);
    let tet_points = BVec::from_iter_in(verts.iter().map(|vid| point(&tets.points, vid.0)), bump);
    let trans_vmat = bumpalo::vec![in bump;
        sub_short::<3, _>(tet_points[1], tet_points[0]),
        sub_short::<3, _>(tet_points[2], tet_points[0]),
        sub_short::<3, _>(tet_points[3], tet_points[0]),
        sub_short::<3, _>(tet_points[2], tet_points[1]),
        sub_short::<3, _>(tet_points[3], tet_points[1]),
        sub_short::<3, _>(tet_points[3], tet_points[2]),
    ];

    let sq_det_vmat = {
        let d = det(&trans_vmat);
        d * d
    };

    let adj_vmat = bumpalo::vec![in bump;
        cross(&trans_vmat[1], &trans_vmat[2]),
        cross(&trans_vmat[2], &trans_vmat[0]),
        cross(&trans_vmat[0], &trans_vmat[1]),
    ];

    let mut interpolant_vec = Vec::with_capacity_in(surfs.len(), bump);
    let mut interpolant_diff_vec = Vec::with_capacity_in(surfs.len(), bump);
    let mut val_diff_vec = Vec::with_capacity_in(surfs.len(), bump);
    for (i, &sid) /*surface id*/ in surfs.iter().enumerate() {
        let tet_vals_grads = BVec::from_iter_in(
            verts.iter().map(|&vid| data.vals_and_grads[sid].get(&vid.0).unwrap()), bump);
        let mut vals = Vec::with_capacity_in(20, bump);

        vals.extend(tet_vals_grads.iter().map(|vals_grads| vals_grads[0]));
        let v0 = tet_vals_grads[0][0];
        let g = &tet_vals_grads[0][1..];
        const S: f64 = 1.0 / 3.0;
        vals.extend([v0 + S * dot(g, &trans_vmat[0]), v0 + S * dot(g, &trans_vmat[1]), v0 + S * dot(g, &trans_vmat[2])]);

        let v1 = tet_vals_grads[1][0];
        let g = &tet_vals_grads[1][1..];
        vals.extend([v1 + S *dot(g, &trans_vmat[3]), v1 + S * dot(g, &trans_vmat[4]), v1 - S * dot(g, &trans_vmat[0])]);

        let v2 = tet_vals_grads[2][0];
        let g = &tet_vals_grads[2][1..];
        vals.extend([v2 + S * dot(g, &trans_vmat[5]), v2 - S * dot(g, &trans_vmat[1]), v2 - S * dot(g, &trans_vmat[3])]);


        let v3 = tet_vals_grads[3][0];
        let g = &tet_vals_grads[3][1..];
        vals.extend([v3 - S * dot(g, &trans_vmat[2]), v3 - S * dot(g, &trans_vmat[4]), v3 - S * dot(g, &trans_vmat[5])]);

        vals.push(((vals[7] + vals[8] + vals[10] + vals[12] + vals[14] + vals[15]) * 1.5 - vals[1] - vals[2] - vals[3]) / 6.0);
        vals.push(((vals[5] + vals[6] + vals[10] + vals[11] + vals[13] + vals[15]) * 1.5 - vals[0] - vals[2] - vals[3]) / 6.0);
        vals.push(((vals[4] + vals[6] + vals[8] + vals[9] + vals[13] + vals[14]) * 1.5 - vals[0] - vals[1] - vals[3]) / 6.0);
        vals.push(((vals[4] + vals[5] + vals[7] + vals[9] + vals[11] + vals[12]) * 1.5 - vals[0] - vals[1] - vals[2]) / 6.0);


        let mut diffs = Vec::with_capacity_in(16, bump);
        for i in 0..16 {
            let c = &C[i];
            diffs.push(v0 * c[0] + v1 * c[1] + v2 * c[2] + v3 * c[3] - vals[i + 4]);
        }

        let val_diff = [v1 - v0, v2 - v0, v3 - v0];
        if test_distance_1(&adj_vmat, val_diff, &diffs, sq_det_vmat, sq_eps) {
            let mut iter = active.iter();
            surfs.retain(|_| *iter.next().unwrap());
            return true;
        } else {
            if *vals.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap() <= 0.0 ||
                *vals.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap() >= 0.0 {
                    active[i] = false;
            }
        }

        interpolant_vec.push(vals);
        interpolant_diff_vec.push(diffs);
        val_diff_vec.push(val_diff);
    }

    let activated = BVec::from_iter_in(
        active
            .iter()
            .enumerate()
            .filter(|(_, &v)| v)
            .map(|(i, _)| i),
        bump,
    );
    for (&i, &j) in activated.iter().tuple_combinations() {
        let points = bumpalo::collections::Vec::from_iter_in(
            interpolant_vec[i]
                .iter()
                .interleave(&interpolant_vec[j])
                .map(|v| *v),
            bump,
        );
        if !contain_dim_2(&points, &[0.0, 0.0], bump) {
            continue;
        }

        let h = bumpalo::vec![in bump; val_diff_vec[i], val_diff_vec[j]];
        let b = [
            interpolant_diff_vec[i].as_slice(),
            interpolant_diff_vec[j].as_slice(),
        ];
        if test_distance_2(&adj_vmat, &h, b, sq_det_vmat, sq_eps, bump) {
            let mut iter = active.iter();
            surfs.retain(|_| *iter.next().unwrap());
            return true;
        }
    }

    for (&i, &j, &k) in activated.iter().tuple_combinations() {
        let points = bumpalo::collections::Vec::from_iter_in(
            interpolant_vec[i]
                .iter()
                .zip(&interpolant_vec[j])
                .zip(&interpolant_vec[k])
                .map(|((a, b), c)| [*a, *b, *c])
                .flatten(),
            bump,
        );
        if !contain_dim_3(&points, &[0.0, 0.0, 0.0], bump) {
            continue;
        }

        let h = bumpalo::vec![in bump; val_diff_vec[i], val_diff_vec[j], val_diff_vec[k]];
        let b = [
            interpolant_diff_vec[i].as_slice(),
            interpolant_diff_vec[j].as_slice(),
            interpolant_diff_vec[k].as_slice(),
        ];
        if test_distance_3(&adj_vmat, &h, b, sq_det_vmat, sq_eps, bump) {
            let mut iter = active.iter();
            surfs.retain(|_| *iter.next().unwrap());
            return true;
        }
    }
    {
        let mut iter = active.into_iter();
        surfs.retain(|_| iter.next().unwrap());
    }
    false
}

fn transpose_adjacent_mat<'b, const N: usize>(
    mat: &[[f64; N]],
    bump: &'b Bump,
) -> bumpalo::collections::Vec<'b, [f64; N]> {
    let mut vec = bumpalo::vec![in bump; [0.0; N]; N];
    if N == 2 {
        vec[0][0] = mat[1][1];
        vec[0][1] = -mat[0][1];
        vec[1][0] = -mat[1][0];
        vec[1][1] = mat[0][0];
    } else if N == 3 {
        cross_in(&mat[1], &mat[2], &mut vec[0]);
        cross_in(&mat[2], &mat[0], &mut vec[1]);
        cross_in(&mat[0], &mat[1], &mut vec[2]);
    } else {
        panic!("not implemented");
    }
    vec
}

fn det<const N: usize>(mat: &[[f64; N]]) -> f64 {
    if N == 2 {
        mat[0][0] * mat[1][1] - mat[0][1] * mat[0][1]
    } else if N == 3 {
        mat[0][0] * mat[1][1] * mat[2][2]
            + mat[0][1] * mat[1][2] * mat[2][0]
            + mat[0][2] * mat[1][0] * mat[2][1]
            - mat[0][2] * mat[1][1] * mat[2][0]
            - mat[0][1] * mat[1][0] * mat[2][2]
            - mat[0][0] * mat[1][2] * mat[2][1]
    } else {
        panic!("not implemented");
    }
}

fn contain_dim_2(points: &[f64], query: &[f64], bump: &Bump) -> bool {
    for ((i, pa), (j, pb)) in points.chunks(2).enumerate().tuple_combinations() {
        let base_ori = double_to_sign(orient2d(pa, pb, query, bump));
        if base_ori == Orientation::Zero {
            continue;
        }

        let separating = || {
            for (k, pc) in points.chunks(2).enumerate() {
                if i == k || j == k {
                    continue;
                }
                let ori = double_to_sign(orient2d(pa, pb, pc, bump));
                if ori == base_ori {
                    return false;
                }
            }
            true
        };
        if separating() {
            return false;
        }
    }
    true
}

fn contain_dim_3(points: &[f64], query: &[f64], bump: &Bump) -> bool {
    for ((i, pa), (j, pb), (k, pc)) in points.chunks(3).enumerate().tuple_combinations() {
        let base_ori = double_to_sign(orient3d(pa, pb, pc, query, bump));
        if base_ori == Orientation::Zero {
            continue;
        }

        let separating = || {
            for (l, pd) in points.chunks(3).enumerate() {
                if i == l || j == l || k == l {
                    continue;
                }
                let ori = double_to_sign(orient3d(pa, pb, pc, pd, bump));
                if ori == base_ori {
                    return false;
                }
            }
            true
        };
        if separating() {
            return false;
        }
    }
    true
}

fn test_distance_1(adj_v: &[[f64; 3]], h: [f64; 3], b: &[f64], sq_det_v: f64, sq_eps: f64) -> bool {
    // w: (M, 3)
    let mut w = [0.0f64; 3];
    for i in 0..3 {
        w[i] = h[0] * adj_v[0][i] + h[1] * adj_v[1][i] + h[2] * adj_v[2][i];
    }
    let w2 = square_norm(&w);
    let max_b = b
        .iter()
        .map(|s| s.abs())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let b2 = max_b * max_b;
    return b2 * sq_det_v > w2 * sq_eps;
}

fn test_distance_2(
    adj_v: &[[f64; 3]],
    h: &[[f64; 3]],
    b: [&[f64]; 2],
    sq_det_v: f64,
    sq_eps: f64,
    bump: &Bump,
) -> bool {
    // w: (M, 3)
    let mut w = bumpalo::vec![in bump; [0.0f64; 3]; 2];
    for i in 0..2 {
        for j in 0..3 {
            w[i][j] = h[i][0] * adj_v[0][j] + h[i][1] * adj_v[1][j] + h[i][2] * adj_v[2][j];
        }
    }
    // u = w * w^T with shape(M, M)
    let mut u = bumpalo::vec![in bump; [0.0; 2]; 2];

    for i in 0..2 {
        for j in 0..2 {
            u[i][j] = w[i][0] * w[j][0] + w[i][1] * w[j][1] + w[i][2] * w[j][2];
        }
    }

    let det_u = det(&u);

    // adj_u: (M, M)
    let trans_adj_u = transpose_adjacent_mat(&u, bump);
    // wu = w^T x adj_u with shape (3, M)
    let mut wu = bumpalo::vec![in bump; [0.0; 2]; 3];
    for i in 0..3 {
        for j in 0..2 {
            for k in 0..2 {
                wu[i][j] += w[k][i] * trans_adj_u[j][k];
            }
        }
    }
    let r2 = (0..b[0].len())
        .map(|l| {
            let mut d = [0.0; 3];
            for i in 0..3 {
                for j in 0..2 {
                    d[i] += wu[i][j] * b[j][l];
                }
            }
            square_norm(&d)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let sq_det_u = det_u * det_u;

    return r2 * sq_det_v > sq_det_u * sq_eps;
}

fn test_distance_3(
    adj_v: &[[f64; 3]],
    h: &[[f64; 3]],
    b: [&[f64]; 3],
    sq_det_v: f64,
    sq_eps: f64,
    bump: &Bump,
) -> bool {
    // w: (M, 3)
    let mut w = bumpalo::vec![in bump; [0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            w[i][j] = h[i][0] * adj_v[0][j] + h[i][1] * adj_v[1][j] + h[i][2] * adj_v[2][j];
        }
    }
    let det_w = det(&w);
    let trans_adj_w = transpose_adjacent_mat(&w, bump);
    let r2 = (0..b[0].len())
        .map(|l| {
            let mut d = [0.0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    d[i] += trans_adj_w[i][j] * b[j][l];
                }
            }
            square_norm(&d)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    return r2 * sq_det_v > det_w * det_w * sq_eps;
}
