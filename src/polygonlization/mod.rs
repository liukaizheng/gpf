mod bsp_complex;
mod conforming_mesh;

use bumpalo::{collections::Vec, Bump};

use crate::{
    predicates::get_exponent,
    triangle::{tetrahedralize, triangulate_polygon_soup},
};

use self::{bsp_complex::BSPComplex, conforming_mesh::Constraints};

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

fn make_mesh_for_triangles<'a, 'b: 'a>(points: &'a [f64], triangles: Vec<'b, usize>) {
    let bump = triangles.bump();
    let mut constraints = Constraints::new(triangles);
    let mut tet_mesh = tetrahedralize(points, bump);
    let tet_triangles = Vec::from_iter_in(
        tet_mesh
            .tets
            .iter()
            .enumerate()
            .filter_map(|(tid, t)| {
                if tet_mesh.is_hull_tet(tid) {
                    None
                } else {
                    let d = &t.data;
                    Some([
                        d[0], d[2], d[1], d[0], d[1], d[3], d[1], d[2], d[3], d[2], d[0], d[3],
                    ])
                }
            })
            .flatten(),
        bump,
    );
    constraints.place_virtual_constraints(&tet_mesh);
    let tet_marks = constraints.insert_constraints(&mut tet_mesh);
    let mut complex = BSPComplex::new(tet_mesh, constraints, tet_marks);
    let mut cid = 0;
    let mut split_bump = Bump::new();
    while cid < complex.n_cells() {
        if complex.splittable(cid) {
            println!("split cell {}", cid);
            split_bump.reset();
            complex.split_cell(cid, &split_bump);
        } else {
            cid += 1;
        }
    }
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
    let (triangles, tri_parents) = triangulate_polygon_soup(&points, &face_edges, axis_data, bump);
    make_mesh_for_triangles(&points, triangles);
}
