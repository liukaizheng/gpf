#![feature(allocator_api)]

use std::path::Path;

use gpf::{boolean3d::{boolean3d, BrepModel}, geometry::{arc::Arc, polyline::Polyline, Crv, Cylinder, Plane, Surf}, math::dot, point};

fn get_angle(plane: &Plane, pt: &[f64]) -> f64 {
    let o = &plane.o;
    let v = [pt[0] - o[0], pt[1] - o[1], pt[2] - o[2]];
    let x = dot(&v, &plane.dx);
    let y = dot(&v, &plane.dy);
    y.atan2(x)
}

fn get_polyline_from_file(name: &str) -> Vec<f64> {
    let file_path = Path::new("tests/data/boolean").join(name);
    let txt = std::fs::read_to_string(file_path).unwrap();
    serde_json::from_str::<Vec<f64>>(&txt).unwrap()
}

fn main() {
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
        1e-6,
    );
}
