use std::collections::HashMap;

use bumpalo::{collections::Vec, vec, Bump};

use crate::math::{dot, sub};
use crate::{predicates, INVALID_IND};

pub fn triangulate<'b>(points: &[f64], segments: &[usize], bump: &'b Bump) -> Vec<'b, usize> {
    let n_points = points.len() >> 1;
    let mut sorted_pt_inds = Vec::from_iter_in(0..n_points, bump);
    sorted_pt_inds.sort_unstable_by(|&i, &j| {
        let i = i << 1;
        let j = j << 1;
        (points[i], points[i + 1])
            .partial_cmp(&(points[j], points[j + 1]))
            .unwrap()
    });
    // resort the array of points to accommodate alternating cuts
    alternate_axes(points, &mut sorted_pt_inds, 0);

    let mut mesh = Mesh {
        points,
        triangles: Vec::new_in(bump),
    };
    let mut hull_left = HEdge::default();
    let mut hull_right = HEdge::default();
    div_conq_recurse(
        &mut mesh,
        &sorted_pt_inds,
        0,
        &mut hull_left,
        &mut hull_right,
        bump,
    );

    let mut ghost = mark_ghost(&mesh.triangles, &mut hull_left, bump);
    form_skeleton(&mut mesh, &ghost, segments, bump);
    let mut visited = vec![in bump; false; mesh.triangles.len()];
    for i in 0..visited.len() {
        if visited[i] || ghost[i] {
            continue;
        }
        visited[i] = true;
        let mut queue: Vec<'_, usize> = vec![in bump; i];
        let mut idx = 0;
        let mut flag = 2;
        while idx < queue.len() {
            let tri = &mesh.triangles[queue[idx]];
            for j in 0..3 {
                if tri.data[j + 3] != INVALID_IND {
                    if flag > 1 {
                        flag = if (tri.data[j + 3] & 1) == 0 { 1 } else { 0 };
                    }
                } else {
                    let he = &tri.nei[j];
                    if visited[he.tri] || ghost[he.tri] {
                        continue;
                    }
                    visited[he.tri] = true;
                    queue.push(he.tri);
                }
            }
            idx += 1;
        }
        if flag == 0 {
            for j in queue {
                ghost[j] = true;
            }
        }
    }
    Vec::from_iter_in(
        mesh.triangles
            .iter()
            .zip(ghost)
            .filter_map(|(tri, is_ghost)| {
                if is_ghost {
                    None
                } else {
                    let data = &tri.data;
                    Some([data[0], data[1], data[2]])
                }
            })
            .flatten(),
        bump,
    )
}

#[inline(always)]
fn point3(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

#[inline]
pub fn triangulate_polygon<'b>(
    points: &[f64],
    segments: &[usize],
    o: &[f64],
    x: &[f64],
    y: &[f64],
    bump: &'b Bump,
) -> Vec<'b, usize> {
    let [new_segments, new_to_ori_map] = unique_indices(segments, bump);

    let points_2d = Vec::from_iter_in(
        new_to_ori_map
            .iter()
            .map(|&idx| {
                let p = point3(points, idx);
                let mut v = vec![in bump; 0.0; 3];
                sub(p, o, &mut v);
                [dot(&v, x), dot(&v, y)]
            })
            .flatten(),
        bump,
    );
    Vec::from_iter_in(
        triangulate(&points_2d, &new_segments, bump)
            .into_iter()
            .map(|idx| new_to_ori_map[idx]),
        bump,
    )
}

#[inline]
pub fn triangulate_polygon_soup<'b>(
    points: &[f64],
    edges: &[Vec<'b, usize>],
    axes: &[f64],
    bump: &'b Bump,
) -> (Vec<'b, usize>, Vec<'b, usize>) {
    let mut triangles = Vec::new_in(bump);
    let mut parents = Vec::new_in(bump);
    for (idx, (segments, axis_data)) in edges.iter().zip(axes.chunks(9)).enumerate() {
        let face_triangles = triangulate_polygon(
            points,
            segments,
            &axis_data[0..3],
            &axis_data[3..6],
            &axis_data[6..9],
            bump,
        );
        parents.resize(parents.len() + face_triangles.len() / 3, idx);
        triangles.extend(face_triangles);
    }
    (triangles, parents)
}

#[derive(Clone, PartialEq, Eq)]
struct HEdge {
    tri: usize,
    ori: usize,
}

impl Default for HEdge {
    fn default() -> Self {
        Self {
            tri: INVALID_IND,
            ori: 0,
        }
    }
}

struct Triangle {
    data: [usize; 6],
    nei: [HEdge; 3],
}

impl Default for Triangle {
    fn default() -> Self {
        Self {
            data: [INVALID_IND; 6],
            nei: [HEdge::default(), HEdge::default(), HEdge::default()],
        }
    }
}

struct Mesh<'a, 'b: 'a> {
    points: &'a [f64],
    triangles: Vec<'b, Triangle>,
}

fn alternate_axes(points: &[f64], indices: &mut [usize], mut axis: usize) {
    let len = indices.len();
    let divider = len >> 1;
    if len <= 3 {
        axis = 0;
    }
    indices.select_nth_unstable_by(divider, |&i, &j| {
        let pi = &points[(i << 1)..];
        let pj = &points[(j << 1)..];
        (pi[axis], pi[1 - axis])
            .partial_cmp(&(pj[axis], pj[1 - axis]))
            .unwrap()
    });
    let (left, right) = indices.split_at_mut(divider);
    if len - divider >= 2 {
        if divider >= 2 {
            alternate_axes(points, left, 1 - axis);
        }
        alternate_axes(points, right, 1 - axis);
    }
}

#[inline(always)]
fn make_triangle(triangles: &mut Vec<'_, Triangle>, he: &mut HEdge) {
    he.tri = triangles.len();
    triangles.push(Triangle::default());
    he.ori = 0;
}

#[inline(always)]
fn set_org(triangles: &mut [Triangle], he: &HEdge, vid: usize) {
    triangles[he.tri].data[(he.ori + 1) % 3] = vid;
}

#[inline(always)]
fn set_dest(triangles: &mut [Triangle], he: &HEdge, vid: usize) {
    triangles[he.tri].data[(he.ori + 2) % 3] = vid;
}

#[inline(always)]
fn set_apex(triangles: &mut [Triangle], he: &HEdge, vid: usize) {
    triangles[he.tri].data[he.ori] = vid;
}

#[inline(always)]
fn bond(triangles: &mut [Triangle], he1: &HEdge, he2: &HEdge) {
    triangles[he1.tri].nei[he1.ori] = he2.clone();
    triangles[he2.tri].nei[he2.ori] = he1.clone();
}

#[inline(always)]
fn org(triangles: &[Triangle], he: &HEdge) -> usize {
    triangles[he.tri].data[(he.ori + 1) % 3]
}

#[inline(always)]
fn dest(triangles: &[Triangle], he: &HEdge) -> usize {
    triangles[he.tri].data[(he.ori + 2) % 3]
}

#[inline(always)]
fn apex(triangles: &[Triangle], he: &HEdge) -> usize {
    triangles[he.tri].data[he.ori]
}

#[inline(always)]
fn point(points: &[f64], idx: usize) -> &[f64] {
    &points[(idx << 1)..]
}

#[inline(always)]
fn prev_self(he: &mut HEdge) {
    he.ori = (he.ori + 2) % 3;
}

#[inline(always)]
fn next_self(he: &mut HEdge) {
    he.ori = (he.ori + 1) % 3;
}

#[inline(always)]
fn sym_self(triangles: &[Triangle], he: &mut HEdge) {
    let nei = &triangles[he.tri].nei[he.ori];
    he.tri = nei.tri;
    he.ori = nei.ori;
}

#[inline(always)]
fn oprev_self(triangles: &[Triangle], he: &mut HEdge) {
    sym_self(triangles, he);
    next_self(he);
}

#[inline(always)]
fn onext_self(triangles: &[Triangle], he: &mut HEdge) {
    prev_self(he);
    sym_self(triangles, he);
}

#[inline(always)]
fn prev(src: &HEdge, dest: &mut HEdge) {
    dest.tri = src.tri;
    dest.ori = (src.ori + 2) % 3;
}

#[inline(always)]
fn next(src: &HEdge, dest: &mut HEdge) {
    dest.tri = src.tri;
    dest.ori = (src.ori + 1) % 3;
}

#[inline(always)]
fn sym(triangles: &[Triangle], src: &HEdge, dest: &mut HEdge) {
    copy(&triangles[src.tri].nei[src.ori], dest);
}

#[inline(always)]
fn oprev(triangles: &[Triangle], src: &HEdge, dest: &mut HEdge) {
    sym(triangles, src, dest);
    next_self(dest);
}

#[inline(always)]
fn onext(triangles: &[Triangle], src: &HEdge, dest: &mut HEdge) {
    prev(src, dest);
    sym_self(triangles, dest);
}

#[inline(always)]
fn copy(src: &HEdge, dest: &mut HEdge) {
    dest.tri = src.tri;
    dest.ori = src.ori;
}

#[inline(always)]
fn counterclockwise(points: &[f64], i: usize, j: usize, k: usize, bump: &Bump) -> f64 {
    predicates::orient2d(point(points, i), point(points, j), point(points, k), bump)
}

#[inline(always)]
fn incircle(points: &[f64], a: usize, b: usize, c: usize, d: usize, bump: &Bump) -> f64 {
    predicates::incircle(
        point(points, a),
        point(points, b),
        point(points, c),
        point(points, d),
        bump,
    )
}

fn merge_hulls(
    m: &mut Mesh,
    axis: usize,
    far_left: &mut HEdge,
    inner_left: &mut HEdge,
    inner_right: &mut HEdge,
    far_right: &mut HEdge,
    bump: &Bump,
) {
    let mut inner_left_dest = dest(&m.triangles, inner_left);
    let mut inner_left_apex = apex(&m.triangles, inner_left);
    let mut inner_right_org = org(&m.triangles, inner_right);
    let mut inner_right_apex = apex(&m.triangles, inner_right);
    let mut far_left_pt: usize;
    let mut far_left_apex: usize;
    let mut far_right_pt: usize;
    let mut far_right_apex: usize;
    // Special treatment for horizontal cuts
    if axis == 1 {
        far_left_pt = org(&m.triangles, far_left);
        far_left_apex = apex(&m.triangles, far_left);
        far_right_pt = dest(&m.triangles, far_right);
        // The pointers to the extreme vertices are shifted to point to the
        // topmost and bottommost vertex of each hull, rather than the
        // leftmost and rightmost vertices.
        while point(m.points, far_left_apex)[1] < point(m.points, far_left_pt)[1] {
            next_self(far_left);
            sym_self(&mut m.triangles, far_left);
            far_left_pt = far_left_apex;
            far_left_apex = apex(&m.triangles, far_left);
        }

        let mut check_edge = HEdge::default();
        sym(&mut m.triangles, inner_left, &mut check_edge);
        let mut check_vertex = apex(&m.triangles, &mut check_edge);
        while point(m.points, check_vertex)[1] > point(m.points, inner_left_dest)[1] {
            next(&mut check_edge, inner_left);
            inner_left_apex = inner_left_dest;
            inner_left_dest = check_vertex;
            sym(&mut m.triangles, inner_left, &mut check_edge);
            check_vertex = apex(&m.triangles, &check_edge);
        }

        while point(m.points, inner_right_apex)[1] < point(m.points, inner_right_org)[1] {
            next_self(inner_right);
            sym_self(&mut m.triangles, inner_right);
            inner_right_org = inner_right_apex;
            inner_right_apex = apex(&m.triangles, inner_right);
        }

        sym(&mut m.triangles, far_right, &mut check_edge);
        check_vertex = apex(&m.triangles, &mut check_edge);
        while point(m.points, check_vertex)[1] > point(m.points, far_right_pt)[1] {
            next(&check_edge, far_right);
            far_right_pt = check_vertex;
            sym(&mut m.triangles, far_right, &mut check_edge);
            check_vertex = apex(&m.triangles, &check_edge);
        }
    }

    // Find a line tangent to and below both hulls
    let mut change_made = true;
    while change_made {
        change_made = false;
        // Make innerleftdest the "bottommost" vertex of the left hull
        if counterclockwise(
            m.points,
            inner_left_dest,
            inner_left_apex,
            inner_right_org,
            bump,
        ) > 0.0
        {
            prev_self(inner_left);
            sym_self(&mut m.triangles, inner_left);
            inner_left_dest = inner_left_apex;
            inner_left_apex = apex(&m.triangles, inner_left);
            change_made = true;
        }
        // Make innerrightorg the "bottommost" vertex of the right hull
        if counterclockwise(
            m.points,
            inner_right_apex,
            inner_right_org,
            inner_left_dest,
            bump,
        ) > 0.0
        {
            next_self(inner_right);
            sym_self(&mut m.triangles, inner_right);
            inner_right_org = inner_right_apex;
            inner_right_apex = apex(&m.triangles, inner_right);
            change_made = true;
        }
    }

    // Find the two candidates to be the next "gear tooth"
    let mut left_cand = HEdge::default();
    let mut right_cand = HEdge::default();
    let mut base_edge = HEdge::default();
    sym(&mut m.triangles, inner_left, &mut left_cand);
    sym(&mut m.triangles, inner_right, &mut right_cand);
    // Create the bottom new bounding triangle
    make_triangle(&mut m.triangles, &mut base_edge);
    // Connect it to the bounding boxes of the left and right triangulations.
    bond(&mut m.triangles, &base_edge, inner_left);
    next_self(&mut base_edge);
    bond(&mut m.triangles, &base_edge, inner_right);
    next_self(&mut base_edge);
    set_org(&mut m.triangles, &base_edge, inner_right_org);
    set_dest(&mut m.triangles, &base_edge, inner_left_dest);

    // Fix the extreme triangles if necessary
    far_left_pt = org(&m.triangles, far_left);
    if inner_left_dest == far_left_pt {
        next(&base_edge, far_left);
    }
    far_right_pt = dest(&m.triangles, far_right);
    if inner_right_org == far_right_pt {
        prev(&base_edge, far_right);
    }

    // The vertices of the current knitting edge
    let mut lower_left = inner_left_dest;
    let mut lower_right = inner_right_org;
    // The candidate vertices for knitting
    let mut upper_left = apex(&m.triangles, &left_cand);
    let mut upper_right = apex(&m.triangles, &right_cand);
    // Walk up the gap between the two triangulations, knitting them together
    loop {
        let left_finished =
            counterclockwise(m.points, upper_left, lower_left, lower_right, bump) <= 0.0;
        let right_finished =
            counterclockwise(m.points, upper_right, lower_left, lower_right, bump) <= 0.0;
        let mut check_edge = HEdge::default();
        let mut next_edge = HEdge::default();
        if left_finished && right_finished {
            // Create the top new bounding triangle
            make_triangle(&mut m.triangles, &mut next_edge);
            set_org(&mut m.triangles, &next_edge, lower_left);
            set_dest(&mut m.triangles, &next_edge, lower_right);
            // Apex is intentionally left INVALID
            // Connect it to the bounding boxes of the two triangulations
            bond(&mut m.triangles, &next_edge, &base_edge);
            next_self(&mut next_edge);
            bond(&mut m.triangles, &next_edge, &right_cand);
            next_self(&mut next_edge);
            bond(&mut m.triangles, &next_edge, &left_cand);

            // Special treatment for horizontal cuts
            if axis == 1 {
                far_left_pt = org(&m.triangles, far_left);
                far_right_pt = dest(&m.triangles, far_right);
                far_right_apex = apex(&m.triangles, far_right);
                sym(&mut m.triangles, far_left, &mut check_edge);
                let mut check_vertex = apex(&m.triangles, &mut check_edge);
                // The pointers to the extremal vertices are restored to the
                // leftmost and rightmost vertices (rather than topmost and
                // bottommost)
                while point(m.points, check_vertex)[0] < point(m.points, far_left_pt)[0] {
                    prev(&check_edge, far_left);
                    far_left_pt = check_vertex;
                    sym(&mut m.triangles, far_left, &mut check_edge);
                    check_vertex = apex(&m.triangles, &check_edge);
                }
                while point(m.points, far_right_apex)[0] > point(m.points, far_right_pt)[0] {
                    prev_self(far_right);
                    sym_self(&mut m.triangles, far_right);
                    far_right_pt = far_right_apex;
                    far_right_apex = apex(&m.triangles, far_right);
                }
            }
            return;
        }

        // Consider eliminating edges from the left triangulation
        if !left_finished {
            // What vertex would be exposed if an edge were deleted
            prev(&left_cand, &mut next_edge);
            sym_self(&mut m.triangles, &mut next_edge);
            let mut next_apex = apex(&m.triangles, &next_edge);
            // If nextapex is INVALID, then no vertex would be exposed; the
            // triangulation would have been eaten right through.
            if next_apex != INVALID_IND {
                // Check whether the edge is Delaunay
                let mut bad_edge = incircle(
                    m.points,
                    lower_left,
                    lower_right,
                    upper_left,
                    next_apex,
                    bump,
                ) > 0.0;
                let mut top_casing = HEdge::default();
                let mut side_casing = HEdge::default();
                let mut outer_casing = HEdge::default();
                while bad_edge {
                    // Eliminate the edge with an edge flip.  As a result, the
                    // left triangulation will have one more boundary triangle.
                    next_self(&mut next_edge);
                    sym(&mut m.triangles, &next_edge, &mut top_casing);
                    next_self(&mut next_edge);
                    sym(&mut m.triangles, &next_edge, &mut side_casing);
                    bond(&mut m.triangles, &next_edge, &mut top_casing);
                    bond(&mut m.triangles, &left_cand, &mut side_casing);
                    next_self(&mut left_cand);
                    sym(&mut m.triangles, &left_cand, &mut outer_casing);
                    prev_self(&mut next_edge);
                    bond(&mut m.triangles, &next_edge, &mut outer_casing);
                    // Correct the vertices to reflect the edge flip
                    set_org(&mut m.triangles, &left_cand, lower_left);
                    set_dest(&mut m.triangles, &left_cand, INVALID_IND);
                    set_apex(&mut m.triangles, &left_cand, next_apex);
                    set_org(&mut m.triangles, &next_edge, INVALID_IND);
                    set_dest(&mut m.triangles, &next_edge, upper_left);
                    set_apex(&mut m.triangles, &next_edge, next_apex);
                    // Consider the newly exposed vertex
                    upper_left = next_apex;
                    // What vertex would be exposed if another edge were deleted?
                    copy(&side_casing, &mut next_edge);
                    next_apex = apex(&m.triangles, &mut next_edge);
                    if next_apex != INVALID_IND {
                        // Check whether the edge is Delaunay
                        bad_edge = incircle(
                            m.points,
                            lower_left,
                            lower_right,
                            upper_left,
                            next_apex,
                            bump,
                        ) > 0.0;
                    } else {
                        // Avoid eating right through the triangulation
                        bad_edge = false;
                    }
                }
            }
        }

        // Consider eliminating edges from the right triangulation
        if !right_finished {
            next(&right_cand, &mut next_edge);
            sym_self(&mut m.triangles, &mut next_edge);
            let mut next_apex = apex(&m.triangles, &next_edge);
            if next_apex != INVALID_IND {
                // Check whether the edge is Delaunay
                let mut bad_edge = incircle(
                    m.points,
                    lower_left,
                    lower_right,
                    upper_right,
                    next_apex,
                    bump,
                ) > 0.0;
                let mut top_casing = HEdge::default();
                let mut side_casing = HEdge::default();
                let mut outer_casing = HEdge::default();
                while bad_edge {
                    prev_self(&mut next_edge);
                    sym(&mut m.triangles, &next_edge, &mut top_casing);
                    prev_self(&mut next_edge);
                    sym(&mut m.triangles, &next_edge, &mut side_casing);
                    bond(&mut m.triangles, &next_edge, &mut top_casing);
                    bond(&mut m.triangles, &right_cand, &mut side_casing);
                    prev_self(&mut right_cand);
                    sym(&mut m.triangles, &right_cand, &mut outer_casing);
                    next_self(&mut next_edge);
                    bond(&mut m.triangles, &next_edge, &mut outer_casing);

                    set_org(&mut m.triangles, &right_cand, INVALID_IND);
                    set_dest(&mut m.triangles, &right_cand, lower_right);
                    set_apex(&mut m.triangles, &right_cand, next_apex);
                    set_org(&mut m.triangles, &next_edge, upper_right);
                    set_dest(&mut m.triangles, &next_edge, INVALID_IND);
                    set_apex(&mut m.triangles, &next_edge, next_apex);

                    upper_right = next_apex;

                    copy(&side_casing, &mut next_edge);
                    next_apex = apex(&m.triangles, &mut next_edge);
                    if next_apex != INVALID_IND {
                        bad_edge = incircle(
                            m.points,
                            lower_left,
                            lower_right,
                            upper_right,
                            next_apex,
                            bump,
                        ) > 0.0;
                    } else {
                        bad_edge = false;
                    }
                }
            }
        }

        if left_finished
            || (!right_finished
                && (incircle(
                    m.points,
                    upper_left,
                    lower_left,
                    lower_right,
                    upper_right,
                    bump,
                ) > 0.0))
        {
            // Knit the triangulations, adding an edge from `lowerleft'
            // to `upperright'
            bond(&mut m.triangles, &base_edge, &right_cand);
            prev(&right_cand, &mut base_edge);
            set_dest(&mut m.triangles, &base_edge, lower_left);
            lower_right = upper_right;
            sym(&mut m.triangles, &base_edge, &mut right_cand);
            upper_right = apex(&m.triangles, &mut right_cand);
        } else {
            // Knit the triangulations, adding an edge from `upperleft'
            // to `lowerright'
            bond(&mut m.triangles, &base_edge, &left_cand);
            next(&left_cand, &mut base_edge);
            set_org(&mut m.triangles, &base_edge, lower_right);
            lower_left = upper_left;
            sym(&mut m.triangles, &base_edge, &mut left_cand);
            upper_left = apex(&m.triangles, &left_cand);
        }
    }
}

fn div_conq_recurse(
    m: &mut Mesh,
    sorted_pt_inds: &[usize],
    axis: usize,
    far_left: &mut HEdge,
    far_right: &mut HEdge,
    bump: &Bump,
) {
    let len = sorted_pt_inds.len();
    if len == 2 {
        make_triangle(&mut m.triangles, far_left);
        set_org(&mut m.triangles, far_left, sorted_pt_inds[0]);
        set_dest(&mut m.triangles, far_left, sorted_pt_inds[1]);

        make_triangle(&mut m.triangles, far_right);
        set_org(&mut m.triangles, far_right, sorted_pt_inds[1]);
        set_dest(&mut m.triangles, far_right, sorted_pt_inds[0]);
        bond(&mut m.triangles, far_left, far_right);

        prev_self(far_left);
        next_self(far_right);
        bond(&mut m.triangles, far_left, far_right);

        prev_self(far_left);
        next_self(far_right);
        bond(&mut m.triangles, far_left, far_right);

        // ensura that the origin of `farleft` is `start`
        prev(far_right, far_left);
    } else if len == 3 {
        let mut midtri = HEdge::default();
        let mut tri1 = HEdge::default();
        let mut tri2 = HEdge::default();
        let mut tri3 = HEdge::default();

        make_triangle(&mut m.triangles, &mut midtri);
        make_triangle(&mut m.triangles, &mut tri1);
        make_triangle(&mut m.triangles, &mut tri2);
        make_triangle(&mut m.triangles, &mut tri3);
        let area = counterclockwise(
            m.points,
            sorted_pt_inds[0],
            sorted_pt_inds[1],
            sorted_pt_inds[2],
            bump,
        );
        if area == 0.0 {
            // Three collinear vertices; the triangulation is two edges
            set_org(&mut m.triangles, &midtri, sorted_pt_inds[0]);
            set_dest(&mut m.triangles, &midtri, sorted_pt_inds[1]);
            set_org(&mut m.triangles, &tri1, sorted_pt_inds[1]);
            set_dest(&mut m.triangles, &tri1, sorted_pt_inds[0]);
            set_org(&mut m.triangles, &tri2, sorted_pt_inds[2]);
            set_dest(&mut m.triangles, &tri2, sorted_pt_inds[1]);
            set_org(&mut m.triangles, &tri3, sorted_pt_inds[1]);
            set_dest(&mut m.triangles, &tri3, sorted_pt_inds[2]);
            // All apices are intentionally left INVALID.
            bond(&mut m.triangles, &midtri, &tri1);
            bond(&mut m.triangles, &tri2, &tri3);
            next_self(&mut midtri);
            prev_self(&mut tri1);
            next_self(&mut tri2);
            prev_self(&mut tri3);
            bond(&mut m.triangles, &midtri, &tri3);
            bond(&mut m.triangles, &tri1, &tri2);
            next_self(&mut midtri);
            prev_self(&mut tri1);
            next_self(&mut tri2);
            prev_self(&mut tri3);
            bond(&mut m.triangles, &midtri, &tri1);
            bond(&mut m.triangles, &tri2, &tri3);

            copy(&tri1, far_left);
            copy(&tri2, far_right);
        } else {
            // The three vertices are not collinear; the triangulation is one
            // triangle, namely `midtri'.
            set_org(&mut m.triangles, &midtri, sorted_pt_inds[0]);
            set_dest(&mut m.triangles, &tri1, sorted_pt_inds[0]);
            set_org(&mut m.triangles, &tri3, sorted_pt_inds[0]);
            if area > 0.0 {
                // The vertices are in counterclockwise order
                set_dest(&mut m.triangles, &midtri, sorted_pt_inds[1]);
                set_org(&mut m.triangles, &tri1, sorted_pt_inds[1]);
                set_dest(&mut m.triangles, &tri2, sorted_pt_inds[1]);
                set_apex(&mut m.triangles, &midtri, sorted_pt_inds[2]);
                set_org(&mut m.triangles, &tri2, sorted_pt_inds[2]);
                set_dest(&mut m.triangles, &tri3, sorted_pt_inds[2]);
            } else {
                // The vertices are in clockwise order
                set_dest(&mut m.triangles, &midtri, sorted_pt_inds[2]);
                set_org(&mut m.triangles, &tri1, sorted_pt_inds[2]);
                set_dest(&mut m.triangles, &tri2, sorted_pt_inds[2]);
                set_apex(&mut m.triangles, &midtri, sorted_pt_inds[1]);
                set_org(&mut m.triangles, &tri2, sorted_pt_inds[1]);
                set_dest(&mut m.triangles, &tri3, sorted_pt_inds[1]);
            }
            // The topology does not depend on how the vertices are ordered
            bond(&mut m.triangles, &midtri, &tri1);
            next_self(&mut midtri);
            bond(&mut m.triangles, &midtri, &tri2);
            next_self(&mut midtri);
            bond(&mut m.triangles, &midtri, &tri3);
            prev_self(&mut tri1);
            next_self(&mut tri2);
            bond(&mut m.triangles, &tri1, &tri2);
            prev_self(&mut tri1);
            prev_self(&mut tri3);
            bond(&mut m.triangles, &tri1, &tri3);
            next_self(&mut tri2);
            prev_self(&mut tri3);
            bond(&mut m.triangles, &tri2, &tri3);
            // Ensure that the origin of `farleft' is start
            copy(&tri1, far_left);
            // Ensure that the destination of `farright' is `start + 2`
            if area > 0.0 {
                copy(&tri2, far_right);
            } else {
                next(far_left, far_right);
            }
        }
    } else {
        let divider = len >> 1;
        let mut innerleft = HEdge::default();
        let mut innerright = HEdge::default();

        let (left, right) = sorted_pt_inds.split_at(divider);

        div_conq_recurse(m, left, 1 - axis, far_left, &mut innerleft, bump);
        div_conq_recurse(m, right, 1 - axis, &mut innerright, far_right, bump);

        merge_hulls(
            m,
            axis,
            far_left,
            &mut innerleft,
            &mut innerright,
            far_right,
            bump,
        );
    }
}

fn mark_ghost<'b>(triangles: &[Triangle], start: &HEdge, bump: &'b Bump) -> Vec<'b, bool> {
    let mut dissolve_edge = start.clone();
    let mut ghost = vec![in bump; false; triangles.len()];
    loop {
        ghost[dissolve_edge.tri] = true;
        next_self(&mut dissolve_edge);
        sym_self(triangles, &mut dissolve_edge);
        if &dissolve_edge == start {
            break;
        }
    }
    ghost
}

fn make_vertex_map<'b>(
    triangles: &[Triangle],
    ghost: &[bool],
    n_points: usize,
    bump: &'b Bump,
) -> Vec<'b, HEdge> {
    let mut vertex_map = vec![in bump; HEdge::default(); n_points];
    for idx in 0..triangles.len() {
        if ghost[idx] {
            continue;
        }
        let tri = &triangles[idx];
        for i in 0..3 {
            let org = tri.data[i];
            if vertex_map[org].tri != INVALID_IND {
                continue;
            }
            let he = &tri.nei[(i + 1) % 3];
            if !ghost[he.tri] {
                copy(he, &mut vertex_map[org]);
            } else {
                vertex_map[org].tri = idx;
                vertex_map[org].ori = (i + 2) % 3;
            }
        }
    }
    return vertex_map;
}

enum Direction {
    WITHIN,
    LEFTCOLLINEAR,
    RIGHTCOLLINEAR,
}

fn find_direction(m: &Mesh, search_tri: &mut HEdge, search_point: usize, bump: &Bump) -> Direction {
    let start_vertex = org(&m.triangles, search_tri);
    let mut right_vertex = dest(&m.triangles, search_tri);
    let mut left_vertex = apex(&m.triangles, search_tri);

    let mut left_ccw = counterclockwise(m.points, search_point, start_vertex, left_vertex, bump);
    let mut left_flag = left_ccw > 0.0;

    let mut right_ccw = counterclockwise(m.points, start_vertex, search_point, right_vertex, bump);
    let mut right_flag = right_ccw > 0.0;
    let mut check_tri = HEdge::default();
    if left_flag && right_flag {
        onext(&m.triangles, search_tri, &mut check_tri);
        if check_tri.tri == INVALID_IND {
            left_flag = false;
        } else {
            right_flag = false;
        }
    }

    while left_flag {
        // Turn left until satisfied
        onext_self(&m.triangles, search_tri);
        left_vertex = apex(&m.triangles, search_tri);
        right_ccw = left_ccw;
        left_ccw = counterclockwise(m.points, search_point, start_vertex, left_vertex, bump);
        left_flag = left_ccw > 0.0;
    }

    while right_flag {
        // Turn right until satisfied
        oprev_self(&m.triangles, search_tri);
        right_vertex = dest(&m.triangles, search_tri);
        left_ccw = right_ccw;
        right_ccw = counterclockwise(m.points, start_vertex, search_point, right_vertex, bump);
        right_flag = right_ccw > 0.0;
    }
    if left_ccw == 0.0 {
        Direction::LEFTCOLLINEAR
    } else if right_ccw == 0.0 {
        Direction::RIGHTCOLLINEAR
    } else {
        Direction::WITHIN
    }
}

#[inline(always)]
fn twin(mark: usize) -> usize {
    if (mark & 1) == 0 {
        mark + 1
    } else {
        mark - 1
    }
}

#[inline(always)]
fn mark(triangles: &[Triangle], he: &HEdge) -> usize {
    triangles[he.tri].data[he.ori + 3]
}
#[inline(always)]
fn set_mark(triangles: &mut [Triangle], he: &HEdge, mark: usize, reverse: bool) {
    triangles[he.tri].data[he.ori + 3] = mark;
    if reverse {
        let mut sym_edge = HEdge::default();
        sym(triangles, he, &mut sym_edge);
        if sym_edge.tri != INVALID_IND {
            triangles[sym_edge.tri].data[sym_edge.ori + 3] = twin(mark);
        }
    }
}

fn scout_segment(
    m: &mut Mesh,
    search_tri: &mut HEdge,
    endpoint2: usize,
    mark: usize,
    bump: &Bump,
) -> bool {
    let collinear = find_direction(m, search_tri, endpoint2, bump);
    let right_vertex = dest(&m.triangles, search_tri);
    let left_vertex = apex(&m.triangles, search_tri);
    if left_vertex == endpoint2 {
        // The segment is already an edge in the mesh.
        prev_self(search_tri);
        set_mark(&mut m.triangles, search_tri, twin(mark), true);
        return true;
    } else if right_vertex == endpoint2 {
        // The segment is already an edge in the mesh.
        set_mark(&mut m.triangles, search_tri, mark, true);
        return true;
    }

    match collinear {
        Direction::WITHIN => false,
        Direction::LEFTCOLLINEAR => {
            prev_self(search_tri);
            scout_segment(m, search_tri, endpoint2, twin(mark), bump)
        }
        Direction::RIGHTCOLLINEAR => {
            next_self(search_tri);
            scout_segment(m, search_tri, endpoint2, mark, bump)
        }
    }
}

fn flip(m: &mut Mesh, flip_edge: &HEdge, vertex_map: &mut [HEdge]) {
    let right_vertex = org(&m.triangles, flip_edge);
    let left_vertex = dest(&m.triangles, flip_edge);
    let bot_vertex = apex(&m.triangles, flip_edge);
    let mut top = HEdge::default();
    sym(&m.triangles, flip_edge, &mut top);
    let far_vertex = apex(&m.triangles, &mut top);

    // Identify the casing of the quadrilateral.
    let mut top_left = HEdge::default();
    prev(&top, &mut top_left);
    let mut topl_casing = HEdge::default();
    sym(&m.triangles, &top_left, &mut topl_casing);
    let mut top_right = HEdge::default();
    next(&top, &mut top_right);
    let mut topr_casing = HEdge::default();
    sym(&m.triangles, &top_right, &mut topr_casing);
    let mut bot_left = HEdge::default();
    next(flip_edge, &mut bot_left);
    let mut botl_casing = HEdge::default();
    sym(&m.triangles, &bot_left, &mut botl_casing);
    let mut bot_right = HEdge::default();
    prev(flip_edge, &mut bot_right);
    let mut botr_casing = HEdge::default();
    sym(&m.triangles, &bot_right, &mut botr_casing);
    // Rotate the quadrilateral one-quarter turn counterclockwise
    bond(&mut m.triangles, &top_left, &botl_casing);
    bond(&mut m.triangles, &bot_left, &botr_casing);
    bond(&mut m.triangles, &bot_right, &topr_casing);
    bond(&mut m.triangles, &top_right, &topl_casing);

    let tl_mark = mark(&m.triangles, &top_left);
    let tr_mark = mark(&m.triangles, &top_right);
    let bl_mark = mark(&m.triangles, &bot_left);
    let br_mark = mark(&m.triangles, &bot_left);

    if &vertex_map[right_vertex] == flip_edge || vertex_map[right_vertex] == top_right {
        copy(&bot_right, &mut vertex_map[right_vertex]);
    }
    if &vertex_map[left_vertex] == &top || vertex_map[left_vertex] == bot_left {
        copy(&top_left, &mut vertex_map[left_vertex]);
    }
    if &vertex_map[bot_vertex] == &bot_right {
        copy(&bot_left, &mut vertex_map[bot_vertex]);
    }
    if &vertex_map[far_vertex] == &top_left {
        copy(&top_right, &mut vertex_map[far_vertex]);
    }

    set_mark(&mut m.triangles, &top_right, tl_mark, false);
    set_mark(&mut m.triangles, &bot_right, tr_mark, false);
    set_mark(&mut m.triangles, &top_left, bl_mark, false);
    set_mark(&mut m.triangles, &bot_left, br_mark, false);

    set_org(&mut m.triangles, flip_edge, far_vertex);
    set_dest(&mut m.triangles, flip_edge, bot_vertex);
    set_apex(&mut m.triangles, flip_edge, right_vertex);
    set_org(&mut m.triangles, &top, bot_vertex);
    set_dest(&mut m.triangles, &top, far_vertex);
    set_apex(&mut m.triangles, &top, left_vertex);
}

fn delaunay_fixup(
    m: &mut Mesh,
    ghost: &[bool],
    fixup_tri: &mut HEdge,
    left_side: bool,
    vertex_map: &mut [HEdge],
    bump: &Bump,
) {
    let mut near_tri = HEdge::default();
    next(fixup_tri, &mut near_tri);
    let mut far_tri = HEdge::default();
    sym(&m.triangles, &near_tri, &mut far_tri);
    // Check if the edge opposite the origin of fixuptri can be flipped
    if ghost[far_tri.tri] {
        return;
    }

    if mark(&m.triangles, &near_tri) != INVALID_IND {
        return;
    }
    let near_vertex = apex(&m.triangles, &near_tri);
    let left_vertex = org(&m.triangles, &near_tri);
    let right_vertex = dest(&m.triangles, &near_tri);
    let far_vertex = apex(&m.triangles, &far_tri);
    // Check whether the previous polygon vertex is a reflex vertex
    if left_side {
        if counterclockwise(&m.points, near_vertex, left_vertex, far_vertex, bump) <= 0.0 {
            // leftvertex is a reflex vertex too, nothing can
            // be done until a convex section is found.
            return;
        }
    } else {
        if counterclockwise(&m.points, far_vertex, right_vertex, near_vertex, bump) <= 0.0 {
            return;
        }
    }

    if counterclockwise(&m.points, right_vertex, left_vertex, far_vertex, bump) > 0.0 {
        if incircle(
            &m.points,
            left_vertex,
            far_vertex,
            right_vertex,
            near_vertex,
            bump,
        ) <= 0.0
        {
            return;
        }
    }

    flip(m, &near_tri, vertex_map);
    prev_self(fixup_tri);
    delaunay_fixup(m, ghost, fixup_tri, left_side, vertex_map, bump);
    delaunay_fixup(m, ghost, &mut far_tri, left_side, vertex_map, bump);
}

fn constrained_edge(
    m: &mut Mesh,
    ghost: &[bool],
    start_tri: &mut HEdge,
    endpoint2: usize,
    mark: usize,
    vertex_map: &mut [HEdge],
    bump: &Bump,
) {
    let endpoint1 = org(&m.triangles, start_tri);
    let mut fixup_tri = HEdge::default();
    next(start_tri, &mut fixup_tri);
    flip(m, &mut fixup_tri, vertex_map);

    let mut collision = false;
    let mut done = false;
    loop {
        let far_vertex = org(&m.triangles, &fixup_tri);
        if far_vertex == endpoint2 {
            let mut fixup_tri2 = HEdge::default();
            oprev(&m.triangles, &fixup_tri, &mut fixup_tri2);
            // Enforce the Delaunay condition around endpoint2
            delaunay_fixup(m, ghost, &mut fixup_tri, false, vertex_map, bump);
            delaunay_fixup(m, ghost, &mut fixup_tri2, true, vertex_map, bump);
            done = true;
        } else {
            let area = counterclockwise(&m.points, endpoint1, endpoint2, far_vertex, bump);
            if area == 0.0 {
                collision = true;
                let mut fixup_tri2 = HEdge::default();
                oprev(&m.triangles, &fixup_tri, &mut fixup_tri2);

                delaunay_fixup(m, ghost, &mut fixup_tri, false, vertex_map, bump);
                delaunay_fixup(m, ghost, &mut fixup_tri2, true, vertex_map, bump);

                done = true;
            } else {
                if area > 0.0 {
                    let mut fixup_tri2 = HEdge::default();
                    oprev(&m.triangles, &fixup_tri, &mut fixup_tri2);
                    delaunay_fixup(m, ghost, &mut fixup_tri2, true, vertex_map, bump);
                    prev_self(&mut fixup_tri);
                } else {
                    delaunay_fixup(m, ghost, &mut fixup_tri, false, vertex_map, bump);
                    oprev_self(&m.triangles, &mut fixup_tri);
                }
                flip(m, &fixup_tri, vertex_map);
            }
        }
        if done {
            break;
        }
    }
    if collision {
        if !scout_segment(m, &mut fixup_tri, endpoint2, mark, bump) {
            constrained_edge(m, ghost, &mut fixup_tri, endpoint2, mark, vertex_map, bump);
        }
    } else {
        let v = dest(&m.triangles, &fixup_tri);
        set_mark(
            &mut m.triangles,
            &fixup_tri,
            if v == endpoint2 { mark } else { twin(mark) },
            true,
        );
    }
}

fn insert_segment(
    m: &mut Mesh,
    vertex_map: &mut [HEdge],
    ghost: &[bool],
    mut start: usize,
    mut end: usize,
    mark: usize,
    bump: &Bump,
) {
    let mut searchtri1 = vertex_map[start].clone();
    if scout_segment(m, &mut searchtri1, end, mark, bump) {
        return;
    }

    // The first endpoint may have changed if a collision with an intervening
    // vertex on the segment occurred.
    start = org(&m.triangles, &searchtri1);

    let mut searchtri2 = vertex_map[end].clone();
    if scout_segment(m, &mut searchtri2, start, twin(mark), bump) {
        return;
    }
    end = org(&m.triangles, &searchtri2);
    constrained_edge(m, ghost, &mut searchtri1, end, mark, vertex_map, bump);
}

fn form_skeleton(m: &mut Mesh, ghost: &[bool], segment: &[usize], bump: &Bump) {
    let mut vertex_map = make_vertex_map(&m.triangles, ghost, m.points.len() >> 1, bump);
    for (i, seg) in segment.chunks(2).enumerate() {
        if seg[0] != seg[1] {
            insert_segment(m, &mut vertex_map, ghost, seg[0], seg[1], i << 1, bump);
        }
    }
}

fn unique_indices<'b>(indices: &[usize], bump: &'b Bump) -> [Vec<'b, usize>; 2] {
    let mut count = 0;
    let mut map = HashMap::new();
    map.reserve(indices.len());
    let mut result = Vec::new_in(bump);
    result.reserve(indices.len());
    for old in indices {
        if let Some(&now) = map.get(old) {
            result.push(now);
        } else {
            map.insert(*old, count);
            result.push(count);
            count += 1;
        }
    }
    let mut new_to_ori_map = vec![in bump; 0; map.len()];
    for (k, v) in map {
        new_to_ori_map[v] = k;
    }
    [result, new_to_ori_map]
}
