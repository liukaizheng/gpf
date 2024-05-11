use gpf::graphcut::{MaxFlow, PushRelabelFifo};

#[test]
fn test_max_flow() {
    let arcs = vec![
        (0, 1, 10),
        (0, 2, 6),
        (1, 3, 4),
        (3, 1, 4),
        (1, 4, 4),
        (4, 1, 4),
        (2, 5, 4),
        (5, 2, 4),
        (3, 4, 2),
        (4, 3, 2),
        (4, 5, 2),
        (5, 4, 2),
        (3, 6, 6),
        (5, 6, 6),
    ];
    let mut graph = PushRelabelFifo::from((arcs, 7));
    graph.find_max_flow();

    println!("haha");
}
