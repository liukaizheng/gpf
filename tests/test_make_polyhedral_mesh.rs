use gpf::polygonlization::{make_polyhedral_mesh, make_mesh_for_triangles};
use serde::Deserialize;

#[allow(non_snake_case)]
#[derive(Deserialize)]
struct PolygonsData {
    pointsData: Vec<f64>,
    faceInShellData: Vec<usize>,
    separators: Vec<usize>,
    edgesData: Vec<usize>,
    axesData: Vec<f64>,
}

fn make_two_dim_arr<T: Clone>(arr: &[T], separators: &[usize]) -> Vec<Vec<T>> {
    Vec::from_iter(
        separators
            .windows(2)
            .into_iter()
            .map(|a| Vec::from_iter(arr[a[0]..a[1]].iter().map(|t| t.clone()))),
    )
}

fn write_obj(points: &[f64], triangles: &[usize], name: &str) {
    let mut txt = "".to_owned();
    for p in points.chunks(3) {
        txt.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
    }
    for tri in triangles.chunks(3) {
        txt.push_str(&format!("f {} {} {}\n", tri[0] + 1, tri[1] + 1, tri[2] + 1));
    }
    std::fs::write(name, txt).unwrap();
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


#[test]
fn two_models() {
    let f = std::fs::File::open("tests/data/solidFix.json").expect("read file");
    let reader = std::io::BufReader::new(f);
    let data = serde_json::from_reader::<_, PolygonsData>(reader).expect("failed to read");
    let edges = make_two_dim_arr(&data.edgesData, &data.separators);
    let (points, triangles) = make_polyhedral_mesh(
        &data.pointsData,
        &data.axesData,
        &data.faceInShellData,
        &edges,
        1e-6,
    );
    write_obj(&points, &triangles, "123.obj");
}

#[test]
fn test_make_polyhedra_mesh() {
    let points = vec![
        0.0, 0.0, 0.0, // 0
        1.0, 0.0, 0.0, // 1
        1.0, 1.0, 0.0, // 2
        0.0, 1.0, 0.0, // 3
        0.0, 0.0, 1.0, // 4
        1.0, 0.0, 1.0, // 5
        1.0, 1.0, 1.0, // 6
        0.0, 1.0, 1.0, // 7
    ];

    let edges = vec![
        vec![0, 3, 3, 2, 2, 1, 1, 0], // bottom
        vec![4, 5, 5, 6, 6, 7, 7, 4], // top
        vec![0, 4, 4, 7, 7, 3, 3, 0], // left
        vec![1, 2, 2, 6, 6, 5, 5, 1], // right
        vec![0, 1, 1, 5, 5, 4, 4, 0], // front
    ];

    let axis = vec![
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, // bottom
        0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // top
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, // left
        1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, // right
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, // front
    ];

    let poly_in_shell = vec![0, 0, 0, 0, 0];
    let (points, triangles) = make_polyhedral_mesh(&points, &axis, &poly_in_shell, &edges, 1e-6);
    write_obj(&points, &triangles, "124.obj");
}

#[test]
fn test_cube_and_sphere() {
    // read cube and sphere
    let (cube_points, cube_tris) = read_obj("tests/data/cube.obj");
    let (sphere_points, sphere_tris) = read_obj("tests/data/sphere.obj");

    // merge cube and sphere into one mesh
    let n_cube_points = cube_points.len() / 3;
    // merge points
    let points = Vec::from_iter(cube_points.into_iter().chain(sphere_points));
    let tri_in_shells = Vec::from_iter(
        // cube is the first shell
        vec![0; cube_tris.len() / 3]
            .into_iter()
            // sphere is second shell
            .chain(vec![1; sphere_tris.len() / 3]),
    );
    // merge triangles
    let triangles = Vec::from_iter(
        cube_tris
            .into_iter()
            .chain(sphere_tris.into_iter().map(|idx| idx + n_cube_points)),
    );
    let (new_points, new_triangles) =
        make_mesh_for_triangles(&points, &triangles, &tri_in_shells);
    write_obj(&new_points, &new_triangles, "125.obj");
}

#[test]
fn test_pig_model() {
    let (points, triangles) = read_obj("tests/data/boarwindmeter.obj");
    let tri_in_shells = vec![0; triangles.len() / 3];
    let (new_points, new_triangles) =
        make_mesh_for_triangles(&points, &triangles, &tri_in_shells);
    write_obj(&new_points, &new_triangles, "126.obj");

}
