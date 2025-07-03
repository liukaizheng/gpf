#![feature(array_chunks)]
#![feature(allocator_api)]
use std::path::Path;

use bumpalo::Bump;
use gpf::{
    boolean3d::{BrepModel, boolean3d},
    geometry::{Crv, Cylinder, Plane, Surf, Surface, arc::Arc, polyline::Polyline},
    math::{cross, dot, normalize, square_norm, sub_short},
    mesh::{EdgeId, ElementId, FaceId, HalfedgeId, ManifoldMesh, Mesh},
    point,
};

fn get_angle(plane: &Plane, pt: &[f64]) -> f64 {
    let o = &plane.o;
    let v = [pt[0] - o[0], pt[1] - o[1], pt[2] - o[2]];
    let x = dot(&v, &plane.dx);
    let y = dot(&v, &plane.dy);
    y.atan2(x)
}

fn read_obj(name: &str) -> (Vec<f64>, Vec<usize>) {
    let (models, _) =
        tobj::load_obj(name, &tobj::LoadOptions::default()).expect("Failed to load obj file");
    let model = &models[0];
    let points = model
        .mesh
        .positions
        .iter()
        .map(|x| *x as f64)
        .collect::<Vec<_>>();
    let triangles = Vec::from_iter(model.mesh.indices.iter().map(|x| *x as usize));
    (points, triangles)
}

fn point_dist_to_srf(srf: &impl Surface, pt: &[f64]) -> f64 {
    srf.eval(pt)[0].abs()
}

#[test]
fn test1() {
    let sq_tol = 1e-10;
    let (points, triangles) = read_obj("data/mesh/shell_22.obj");
    let bump = Bump::new();
    let mesh = ManifoldMesh::new(triangles.array_chunks::<3>(), &bump);

    let start_pt = [0.955037, -0.5, -0.204065];
    let end_pt = [0.563725, -0.5, -0.5];
    // let start_vid = *mesh
    //     .vertices()
    //     .find(|v| {
    //         let p = point::<3>(&points, v.index());
    //         square_norm(&sub_short::<3, _>(p, &start_pt)) < sq_tol
    //     })
    //     .unwrap();
    let start_vid = 2879.into();

    let mut polyline = Vec::from_iter(start_pt);
    let mut curr_vid = start_vid;
    let mut prev_eid = EdgeId::default();
    let srf1 = Cylinder::new(-0.5, 0.0, 0.0, 0.9, 0.3, 0.1f64.sqrt(), 1.1);
    let srf2 = Plane::new(0.0, -0.5, 0.0, 0.0, 1.0, 0.0);
    let compute_normal = |fid: FaceId| -> [f64; 3] {
        let mut he = mesh.face(fid).halfedge();
        let v1 = *he.to();
        he = he.next();
        let v2 = *he.to();
        he = he.next();
        let v3 = *he.to();
        let p1 = point::<3>(&points, *v1);
        let p2 = point::<3>(&points, *v2);
        let p3 = point::<3>(&points, *v3);
        let mut n = cross(&sub_short::<3, _>(p2, p1), &sub_short::<3, _>(p3, p1));
        normalize::<3>(&mut n);
        return n;
    };
    loop {
        let mut curr_hid = HalfedgeId::default();
        let mut finished = false;
        for he in mesh.vertex(curr_vid).outgoing_halfedges() {
            if *he.edge() == prev_eid {
                continue;
            }
            let next_vid = *he.to();
            let pt = point::<3>(&points, *next_vid);
            let d1 = point_dist_to_srf(&srf1, pt);
            let d2 = point_dist_to_srf(&srf2, pt);
            if d1 < 1e-3 && d2 < 1e-3 {
                let n1 = compute_normal(*he.face());
                let n2 = compute_normal(*he.twin().face());
                let dv = dot(&n1, &n2);
                if dv < 0.8 {
                    polyline.extend_from_slice(pt);
                    curr_hid = *he;
                    curr_vid = next_vid;
                    let dist = square_norm(&sub_short::<3, _>(pt, &end_pt));
                    if dist < sq_tol {
                        finished = true;
                    }
                    break;
                }
            }
        }
        if !curr_hid.valid() || finished {
            break;
        } else {
            prev_eid = mesh.he_edge(curr_hid);
        }
    }
    let txt = serde_json::to_string(&polyline).unwrap();
    std::fs::write("polyline.json", txt).unwrap();
    println!("the polyline is {:?}", polyline);
}

fn get_polyline_from_file(name: &str) -> Vec<f64> {
    let file_path = Path::new("tests/data/boolean").join(name);
    let txt = std::fs::read_to_string(file_path).unwrap();
    serde_json::from_str::<Vec<f64>>(&txt).unwrap()
}

#[test]
fn test_boolean1() {
    let model1 = {
        let face_surfaces = vec![
            (
                Surf::Plane(Plane::new(0.0, -0.5, 0.0, 0.0, -1.0, 0.0)),
                false,
            ),
            (
                Surf::Cylinder(Cylinder::new(-0.5, 0.0, 0.0, 0.9, 0.3, 0.1f64.sqrt(), 1.1)),
                false,
            ),
            (Surf::Plane(Plane::new(0.0, 0.5, 0.0, 0.0, 1.0, 0.0)), false),
            (
                Surf::Plane(Plane::new(-0.5, 0.0, 0.0, -1.0, 0.0, 0.0)),
                false,
            ),
            (
                Surf::Plane(Plane::new(0.0, 0.0, -0.5, 0.0, 0.0, -1.0)),
                false,
            ),
            (
                Surf::Cylinder(Cylinder::new(0.0, 0.0, -0.5, 0.0, 1.0, 0.0, 1.0)),
                false,
            ),
        ];
        #[rustfmt::skip]
        let points = vec![
            -0.5, -0.5, -0.5,
            -0.5, -0.5, 0.365815,
            0.9550368189811707, -0.5, -0.20406539738178253,
            0.5637247562408447, -0.5, -0.5,
            1.0, -0.15287435054779053, -0.5,
            -0.5, 0.5, -0.5,
            -0.5, 0.5, 0.365815,
            1.0, 0.5, -0.5
        ];
        let edge_crv_21 = {
            let plane = Plane::from_x_y([0.0, -0.5, -0.5], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]);
            let start_angle = get_angle(&plane, point::<3>(&points, 2));
            let end_angle = get_angle(&plane, point::<3>(&points, 1));
            Crv::Arc(Arc::new(plane, 1.0, start_angle, end_angle))
        };
        let edge_curv_67 = {
            let plane = Plane::from_x_y([0.0, 0.5, -0.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
            let start_angle = get_angle(&plane, point::<3>(&points, 6));
            let end_angle = get_angle(&plane, point::<3>(&points, 7));
            Crv::Arc(Arc::new(plane, 1.0, start_angle, end_angle))
        };
        let edge_curv_34 = Crv::Polyline(Polyline::new(get_polyline_from_file("polyline_34.json")));
        let edge_curv_42 = Crv::Polyline(Polyline::new(get_polyline_from_file("polyline_42.json")));
        let edge_curv_23 = Crv::Polyline(Polyline::new(get_polyline_from_file("polyline_23.json")));
        let edge_curves = vec![
            ([2, 1], edge_crv_21),
            ([6, 7], edge_curv_67),
            ([3, 4], edge_curv_34),
            ([4, 2], edge_curv_42),
            ([2, 3], edge_curv_23),
        ];
        #[rustfmt::skip]
        let loops = vec![
            vec![0, 3, 2, 1],
            vec![2, 3, 4],
            vec![5, 6, 7],
            vec![0, 1, 6, 5],
            vec![5, 7, 4, 3, 0],
            vec![1, 2, 4, 7, 6],
        ];
        let face_loops = vec![vec![0], vec![1], vec![2], vec![3], vec![4], vec![5]];
        BrepModel::new_in(
            loops,
            face_loops,
            points,
            face_surfaces,
            edge_curves,
            std::alloc::Global,
        )
    };

    let model2 = {
        #[rustfmt::skip]
        let points = vec![
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0
        ];
        #[rustfmt::skip]
        let loops = [
            vec![0, 4, 7 ,3],
            vec![1, 2, 6, 5],
            vec![0, 3, 2, 1],
            vec![4, 5, 6, 7],
            vec![0, 1, 5, 4],
            vec![3, 7, 6, 2],
        ];
        let face_loops = vec![vec![0], vec![1], vec![2], vec![3], vec![4], vec![5]];
        let face_surfaces = vec![
            (
                Surf::Plane(Plane::new(0.0, 0.0, 0.0, -1.0, 0.0, 0.0)),
                false,
            ),
            (Surf::Plane(Plane::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)), false),
            (
                Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0)),
                false,
            ),
            (Surf::Plane(Plane::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)), false),
            (
                Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, -1.0, 0.0)),
                false,
            ),
            (Surf::Plane(Plane::new(0.0, 1.0, 0.0, 0.0, 1.0, 0.0)), false),
        ];

        BrepModel::new_in(
            loops,
            face_loops,
            points,
            face_surfaces,
            Vec::new(),
            std::alloc::Global,
        )
    };

    boolean3d(
        vec![model1, model2],
        |is_kept_arr| is_kept_arr[0] && !is_kept_arr[1],
        0.004,
    );
}
