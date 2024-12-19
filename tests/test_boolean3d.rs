use gpf::{
    boolean3d::{boolean3d, BrepModel},
    geometry::{BBox, Cylinder, Plane, Surf},
};

#[test]
fn test_boolean1() {
    let surfaces = vec![
        Surf::Plane(Plane::new(-0.5, 0.0, 0.0, -1.0, 0.0, 0.0)),
        Surf::Cylinder(Cylinder::new(-0.5, 0.0, 0.0, 0.9, 0.3, 0.1f64.sqrt(), 1.1)),
        Surf::Plane(Plane::new(0.0, 0.0, -0.5, 0.0, 0.0, -1.0)),
        Surf::Cylinder(Cylinder::new(0.0, 0.0, -0.5, 0.0, 1.0, 0.0, 1.0)),
        Surf::Plane(Plane::new(0.0, -0.5, 0.0, 0.0, -1.0, 0.0)),
        Surf::Plane(Plane::new(0.0, 0.5, 0.0, 0.0, 1.0, 0.0)),
        Surf::Plane(Plane::new(0.0, 0.0, 0.0, -1.0, 0.0, 0.0)),
        Surf::Plane(Plane::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0)),
        Surf::Plane(Plane::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)),
        Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, -1.0, 0.0)),
        Surf::Plane(Plane::new(0.0, 1.0, 0.0, 0.0, 1.0, 0.0)),
    ];
    let model1 = {
        #[rustfmt::skip]
        let points = vec![
            -0.5, -0.5, -0.5,
            -0.5, -0.5, 0.365815,
            0.955037, -0.5, -0.204065,
            0.563725, -0.5, -0.5,
            1.0, -0.152874, -0.5,
            -0.5, 0.5, -0.5,
            -0.5, 0.5, 0.365815,
            1.0, 0.5, -0.5
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
        let face_surfaces = vec![8, 2, 10, 0, 4, 6];
        BrepModel::new(
            loops,
            face_loops,
            face_surfaces,
            points,
            BBox::new(-0.5, -0.5, -0.5, 1.0, 0.5, 0.365815),
        )
    };

    let model2 = {
        #[rustfmt::skip]
        let points = [
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

        let face_surfaces = vec![12, 14, 16, 18, 20, 22];
        BrepModel::new(
            loops,
            face_loops,
            face_surfaces,
            points,
            BBox::new(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        )
    };

    boolean3d(
        vec![model1, model2],
        surfaces,
        |is_kept_arr| {
            is_kept_arr
                .iter()
                .map(|&e| e)
                .reduce(|res, e| (res | e))
                .unwrap()
        },
        0.1,
    );
}
