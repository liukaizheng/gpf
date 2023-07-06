pub mod bsp_complex;
mod conforming_mesh;

use bumpalo::Bump;

use crate::{
    mesh::Mesh,
    predicates::get_exponent,
    triangle::{tetrahedralize, triangulate_polygon_soup},
};

use self::{bsp_complex::BSPComplex, conforming_mesh::Constraints};

fn point(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

pub fn remove_duplicates(points: &[f64], epsilon: f64) -> (Vec<f64>, Vec<usize>) {
    if epsilon == 0.0 {
        return (
            Vec::from_iter(points.iter().map(|&x| x)),
            Vec::from_iter(0..points.len()),
        );
    }
    let e = get_exponent(epsilon) - 1;
    let base = 2.0f64.powi(e);
    let approx_points = Vec::from_iter(points.iter().map(|x| (x / base).round() * base));
    let n_points = points.len() / 3;
    let mut indices = Vec::from_iter(0..n_points);
    indices.sort_unstable_by(|&i, &j| point(points, i).partial_cmp(point(points, j)).unwrap());
    let mut pmap = vec![0usize; n_points];
    let mut out_pt_indices = vec![indices[0]];
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
        Vec::from_iter(
            out_pt_indices
                .into_iter()
                .map(|idx| {
                    let p = point(&approx_points, idx);
                    [p[0], p[1], p[2]]
                })
                .flatten(),
        ),
        pmap,
    )
}

fn make_mesh_for_triangles(points: &[f64], triangles: Vec<usize>, tri_in_shells: &[usize]) {
    let mut constraints = Constraints::new(triangles);
    let mut tet_mesh = tetrahedralize(points);
    constraints.place_virtual_constraints(&tet_mesh);
    let tet_marks = constraints.insert_constraints(&mut tet_mesh);
    let mut complex = BSPComplex::new(tet_mesh, &constraints, tet_marks);
    {
        let mut cid = 0;
        let mut bump = Bump::new();
        while cid < complex.n_cells() {
            if complex.splittable(cid) {
                bump.reset();
                complex.split_cell(cid, &bump);
            } else {
                cid += 1;
            }
        }
        println!(
            "split cell {}, n_vertices: {}, n_halfedges {}, n_edges {}, n_faces {}",
            cid,
            complex.mesh.n_vertices(),
            complex.mesh.n_halfedges(),
            complex.mesh.n_edges(),
            complex.mesh.n_faces()
        );

        println!("ori is {:?}", complex.ori_duration);
        println!("split is {:?}", complex.split_duration);
    }

    complex.complex_partition(&tri_in_shells);
}

pub fn make_polyhedra_mesh(
    point_data: &[f64],
    axis_data: &[f64],
    face_in_shell_data: &[usize],
    face_edge_data: &[Vec<usize>],
    epsilon: f64,
) {
    let (points, pmap) = remove_duplicates(point_data, epsilon);
    let face_edges = Vec::from_iter(face_edge_data.iter().map(|arr| {
        Vec::from_iter(
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
        )
    }));
    let (triangles, tri_parents) = triangulate_polygon_soup(&points, &face_edges, axis_data);
    make_mesh_for_triangles(
        &points,
        triangles,
        &Vec::from_iter(
            tri_parents
                .into_iter()
                .map(|parent| face_in_shell_data[parent]),
        ),
    );
}
