mod conforming_mesh;

use bumpalo::{collections::Vec, vec, Bump};

use crate::predicates::get_exponent;

fn point(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

pub fn remove_duplicates<'b>(
    points: &[f64],
    bump: &'b Bump,
    epsilon: f64,
) -> (Vec<'b, f64>, Vec<'b, usize>) {
    if epsilon == 0.0 {
        return (
            Vec::from_iter_in(points.iter().map(|&x| x), bump),
            Vec::from_iter_in(0..points.len(), bump),
        );
    }
    let e = get_exponent(epsilon) - 1;
    let base = 2.0f64.powi(e);
    let approx_points = Vec::from_iter_in(points.iter().map(|x| (x / base).round() * base), bump);
    let n_points = points.len() / 3;
    let mut indices = Vec::from_iter_in(0..n_points, bump);
    indices.sort_unstable_by(|&i, &j| point(points, i).partial_cmp(point(points, j)).unwrap());
    let mut pmap = bumpalo::vec![in bump; 0usize; n_points];
    let mut out_pt_indices = bumpalo::vec![in bump; indices[0]];
    for k in 1..n_points {
        let (i, j) = (indices[k], indices[k - 1]);
        if point(points, i) == point(points, j) {
            pmap[i] = out_pt_indices.len() - 1;
        } else {
            pmap[i] = out_pt_indices.len();
            out_pt_indices.push(i);
        }
    }
    (
        Vec::from_iter_in(
            out_pt_indices
                .into_iter()
                .map(|idx| {
                    let p = point(&approx_points, idx);
                    [p[0], p[1], p[2]]
                })
                .flatten(),
            bump,
        ),
        pmap,
    )
}

pub fn make_polyhedra_mesh<'b>(
    point_data: &[f64],
    axis_data: &[f64],
    face_in_shell_data: &[usize],
    face_edge_data: &[Vec<'b, usize>],
    bump: &'b Bump,
    epsilon: f64,
) {
    let (points, pmap) = remove_duplicates(point_data, bump, epsilon);
    let face_edges = Vec::from_iter_in(
        face_edge_data.iter().map(|arr| {
            Vec::from_iter_in(
                arr.chunks(2)
                    .filter_map(|pair| {
                        let a = pmap[pair[0]];
                        let b = pmap[pair[1]];
                        if a == b {
                            None
                        } else {
                            Some([a, b])
                        }
                    })
                    .flatten(),
                bump,
            )
        }),
        bump,
    );
}
