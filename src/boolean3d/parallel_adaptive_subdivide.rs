use std::{default, f64::consts::{SQRT_2, SQRT_3}};

use hashbrown::HashSet;
use itertools::Itertools;
use tinyvec::TinyVec;

use crate::{
    boolean3d::adaptive_subdivide::test_distance_1,
    geometry::{BBox, Surf, Surface},
    math::{cross, dot, interpolate, square_norm, sub_short},
};

use super::adaptive_subdivide::{
    SurfaceData, contain_zero_2, contain_zero_3, det, test_distance_2, test_distance_3,
};

#[derive(Default)]
struct SurfaceEvaluation {
    sid: usize,
    evaluation: [[f64; 4]; 4],
}

#[derive(Default)]
struct Tet {
    points: [[f64; 3]; 4],
    surface_evaluations: TinyVec<[SurfaceEvaluation; 3]>,
    opposite_faces: [TinyVec<[(Box<Tet>, u8); 2]>; 4],
    edge_square_lengths: [f64; 6],
    sub_tets: TinyVec<[(Box<Tet>, u8); 2]>,
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
impl Tet {
    fn subdivide(&mut self, srf_datum: &[SurfaceData], sq_eps: f64) {
        if !self.subdividable(srf_datum, sq_eps) {
            return;
        }

        const T: [[usize; 4]; 6] = [
            [0, 1, 2, 3],
            [0, 2, 3, 1],
            [0, 3, 1, 2],
            [1, 2, 0, 3],
            [1, 3, 2, 0],
            [2, 3, 0, 1],
        ];
        let largest_edge_idx = self
            .edge_square_lengths
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        let table = &T[largest_edge_idx];
        let new_pt = interpolate::<3>(&self.points[table[0]], &self.points[table[1]], 0.5);
        let half_square_length = self.edge_square_lengths[largest_edge_idx] / 4.0;

        let len1 = sub_short::<3, _>(&new_pt, &self.points[table[1]]);
        let len2 = sub_short::<3, _>(&new_pt, &self.points[table[2]]);

        let create_tet = |is_first| {
            let points = [[f64; 3]; 4]::default();
        };
    }
    fn subdividable(&mut self, srf_datum: &[SurfaceData], sq_eps: f64) -> bool {
        let tet_box = BBox::from_iter(&self.points);

        let mut contain_some_srf = false;
        self.surface_evaluations.retain(|&eval| {
            let srf_data = &srf_datum[eval.sid];
            if tet_box.contains(&srf_data.bbox) {
                contain_some_srf = true;
                return true;
            }
            if !tet_box.intersects(&srf_data.bbox) {
                return false;
            }
            if srf_data.sub_bboxes.len() > 1 {
                if srf_data
                    .sub_bboxes
                    .iter()
                    .any(|bbox| tet_box.contains(bbox))
                {
                    contain_some_srf = true;
                    return true;
                }
                for bbox in &srf_data.sub_bboxes {
                    if tet_box.intersects(bbox) {
                        return true;
                    }
                }
                false
            } else {
                true
            }
        });

        if contain_some_srf {
            return true;
        }

        if self.surface_evaluations.len() < 1 {
            return false;
        }
        let tet_points = &self.points;

        let trans_vmat = [
            sub_short::<3, _>(&tet_points[1], &tet_points[0]),
            sub_short::<3, _>(&tet_points[2], &tet_points[0]),
            sub_short::<3, _>(&tet_points[3], &tet_points[0]),
            sub_short::<3, _>(&tet_points[2], &tet_points[1]),
            sub_short::<3, _>(&tet_points[3], &tet_points[1]),
            sub_short::<3, _>(&tet_points[3], &tet_points[2]),
        ];

        let sq_det_vmat = {
            let d = det(&trans_vmat);
            d * d
        };

        let adj_vmat = [
            cross(&trans_vmat[1], &trans_vmat[2]),
            cross(&trans_vmat[2], &trans_vmat[0]),
            cross(&trans_vmat[0], &trans_vmat[1]),
        ];

        let n_surfaces = self.surface_evaluations.len();

        let mut interpolant_vec = Vec::with_capacity(n_surfaces);
        let mut interpolant_diff_vec = Vec::with_capacity(n_surfaces);
        let mut val_diff_vec = Vec::with_capacity(n_surfaces);
        for eval in self.surface_evaluations.iter() {
            let tet_vals_grads = &eval.evaluation;

            let mut vals = Vec::with_capacity(20);
            vals.extend(tet_vals_grads.iter().map(|vals_grads| vals_grads[0]));

            let v0 = tet_vals_grads[0][0];
            let g = &tet_vals_grads[0][1..];
            const S: f64 = 1.0 / 3.0;
            vals.extend([
                v0 + S * dot(g, &trans_vmat[0]),
                v0 + S * dot(g, &trans_vmat[1]),
                v0 + S * dot(g, &trans_vmat[2]),
            ]);

            let v1 = tet_vals_grads[1][0];
            let g = &tet_vals_grads[1][1..];
            vals.extend([
                v1 + S * dot(g, &trans_vmat[3]),
                v1 + S * dot(g, &trans_vmat[4]),
                v1 - S * dot(g, &trans_vmat[0]),
            ]);

            let v2 = tet_vals_grads[2][0];
            let g = &tet_vals_grads[2][1..];
            vals.extend([
                v2 + S * dot(g, &trans_vmat[5]),
                v2 - S * dot(g, &trans_vmat[1]),
                v2 - S * dot(g, &trans_vmat[3]),
            ]);

            let v3 = tet_vals_grads[3][0];
            let g = &tet_vals_grads[3][1..];
            vals.extend([
                v3 - S * dot(g, &trans_vmat[2]),
                v3 - S * dot(g, &trans_vmat[4]),
                v3 - S * dot(g, &trans_vmat[5]),
            ]);

            vals.push(
                ((vals[7] + vals[8] + vals[10] + vals[12] + vals[14] + vals[15]) * 1.5
                    - vals[1]
                    - vals[2]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[5] + vals[6] + vals[10] + vals[11] + vals[13] + vals[15]) * 1.5
                    - vals[0]
                    - vals[2]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[4] + vals[6] + vals[8] + vals[9] + vals[13] + vals[14]) * 1.5
                    - vals[0]
                    - vals[1]
                    - vals[3])
                    / 6.0,
            );
            vals.push(
                ((vals[4] + vals[5] + vals[7] + vals[9] + vals[11] + vals[12]) * 1.5
                    - vals[0]
                    - vals[1]
                    - vals[2])
                    / 6.0,
            );

            let mut diffs = Vec::with_capacity(16);
            for i in 0..16 {
                let c = &C[i];
                diffs.push(v0 * c[0] + v1 * c[1] + v2 * c[2] + v3 * c[3] - vals[i + 4]);
            }

            let val_diff = [v1 - v0, v2 - v0, v3 - v0];
            let is_active = *vals
                .iter()
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap()
                > 0.0
                && *vals
                    .iter()
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap()
                    < 0.0;
            if is_active {
                if test_distance_1(&adj_vmat, val_diff, &diffs, sq_det_vmat, sq_eps) {
                    return true;
                }
                interpolant_vec.push(vals);
                interpolant_diff_vec.push(diffs);
                val_diff_vec.push(val_diff);
            }
        }

        if interpolant_vec.len() < 2 {
            return false;
        }

        let mut pair_set = HashSet::new();
        for ((i, v1), (j, v2)) in interpolant_vec.iter().enumerate().tuple_windows() {
            let mut points = Vec::with_capacity(42);
            points.extend(v1.iter().interleave(v2).map(|v| *v));

            if !contain_zero_2(points, std::alloc::Global) {
                continue;
            }

            pair_set.insert([i, j]);

            let h = [val_diff_vec[i], val_diff_vec[j]];
            let b = [
                interpolant_diff_vec[i].as_slice(),
                interpolant_diff_vec[j].as_slice(),
            ];
            if test_distance_2(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
                return true;
            }
        }

        for ((i, v1), (j, v2), (k, v3)) in interpolant_vec.iter().enumerate().tuple_windows() {
            if !pair_set.contains(&[i, j])
                || !pair_set.contains(&[i, k])
                || !pair_set.contains(&[j, k])
            {
                continue;
            }
            let mut points = Vec::with_capacity(63);
            points.extend(
                v1.iter()
                    .zip(v2)
                    .zip(v3)
                    .map(|((a, b), c)| [*a, *b, *c])
                    .flatten(),
            );

            if !contain_zero_3(points, std::alloc::Global) {
                continue;
            }

            let h = [val_diff_vec[i], val_diff_vec[j], val_diff_vec[k]];
            let b = [
                interpolant_diff_vec[i].as_slice(),
                interpolant_diff_vec[j].as_slice(),
                interpolant_diff_vec[k].as_slice(),
            ];
            if test_distance_3(&adj_vmat, &h, b, sq_det_vmat, sq_eps) {
                return true;
            }
        }
        return false;
    }
}

pub(super) fn build_tet_from_box(bbox: BBox, surfaces: &[Surf]) -> Tet {
    const TET: [[f64; 3]; 4] = [
        [
            1.0 + SQRT_3 / SQRT_2 * 0.5 + 1.0 / SQRT_3,
            -1.0 / (2.0 * SQRT_2),
            0.0,
        ],
        [0.5, 1.0 + (SQRT_2 + SQRT_3) / 2.0, 0.0],
        [
            -1.0 / SQRT_3 - SQRT_3 / SQRT_2 * 0.5,
            -1.0 / (2.0 * SQRT_2),
            0.0,
        ],
        [
            0.5,
            1.0 / 3.0 + SQRT_3 / 6.0,
            1.0 + 2.0 / 3.0 * SQRT_2 + SQRT_2 / SQRT_3,
        ],
    ];
    let len_vec = sub_short::<3, _>(&bbox.max, &bbox.min);
    let tet_points = TET.map(|p| {
        [
            p[0] * len_vec[0] + bbox.min[0],
            p[1] * len_vec[1] + bbox.min[1],
            p[2] * len_vec[2] + bbox.min[2],
        ]
    });

    let edge_square_lengths = [
        square_norm(&sub_short::<3, _>(&tet_points[0], &tet_points[1])),
        square_norm(&sub_short::<3, _>(&tet_points[0], &tet_points[2])),
        square_norm(&sub_short::<3, _>(&tet_points[0], &tet_points[3])),
        square_norm(&sub_short::<3, _>(&tet_points[1], &tet_points[2])),
        square_norm(&sub_short::<3, _>(&tet_points[1], &tet_points[3])),
        square_norm(&sub_short::<3, _>(&tet_points[2], &tet_points[3])),
    ];

    let surface_evaluations =
        TinyVec::from_iter(
            surfaces
                .iter()
                .enumerate()
                .map(|(sid, srf)| SurfaceEvaluation {
                    sid,
                    evaluation: tet_points.map(|p| srf.eval(&p)),
                }),
        );
    Tet {
        points: tet_points,
        surface_evaluations,
        edge_square_lengths,
        opposite_faces: [
            TinyVec::new(),
            TinyVec::new(),
            TinyVec::new(),
            TinyVec::new(),
        ],
        sub_tets: TinyVec::new(),
    }
}

fn write_tet(points_arr: &[[f64; 3]; 4]) {
    let mut file = std::fs::File::create("tet.obj").unwrap();
    use std::io::Write;
    for points in points_arr {
        writeln!(&mut file, "v {} {} {}", points[0], points[1], points[2]).unwrap();
    }

    writeln!(&mut file, "f 1 2 3").unwrap();
    writeln!(&mut file, "f 1 3 4").unwrap();
    writeln!(&mut file, "f 3 2 4").unwrap();
    writeln!(&mut file, "f 2 1 4").unwrap();
}
