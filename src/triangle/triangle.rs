use bumpalo::{collections::Vec, Bump};

use crate::{predicates, INVALID_IND};

pub fn triangulate<'b>(points: &[f64], segments: &[usize], bump: &'b Bump) -> usize {
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

    let mesh = bump.alloc(Mesh {
        points,
        triangles: Vec::new_in(bump),
    });
    let hull_left = bump.alloc(HEdge::default());
    let hull_right = bump.alloc(HEdge::default());
    div_conq_recurse(mesh, &sorted_pt_inds, 0, hull_left, hull_right, bump);
    mesh.triangles.len()
}

#[derive(Clone)]
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
    unsafe {
        *triangles
            .get_unchecked_mut(he.tri)
            .data
            .get_unchecked_mut((he.ori + 1) % 3) = vid;
    }
}

#[inline(always)]
fn set_dest(triangles: &mut [Triangle], he: &HEdge, vid: usize) {
    unsafe {
        *triangles
            .get_unchecked_mut(he.tri)
            .data
            .get_unchecked_mut((he.ori + 2) % 3) = vid;
    }
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
    let nei = &triangles[src.tri].nei[src.ori];
    dest.tri = nei.tri;
    dest.ori = nei.ori;
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
        // The pointers to the extremal vertices are shifted to point to the
        // topmost and bottommost vertex of each hull, rather than the
        // leftmost and rightmost vertices.
        while point(m.points, far_left_apex)[1] < point(m.points, far_left_pt)[1] {
            next_self(far_left);
            sym_self(&mut m.triangles, far_left);
            far_left_pt = far_left_apex;
            far_left_apex = apex(&m.triangles, far_left);
        }

        let check_edge = bump.alloc(HEdge::default());
        sym(&mut m.triangles, inner_left, check_edge);
        let mut check_vertex = apex(&m.triangles, check_edge);
        while point(m.points, check_vertex)[1] > point(m.points, inner_left_dest)[1] {
            next(check_edge, inner_left);
            inner_left_apex = inner_left_dest;
            inner_left_dest = check_vertex;
            sym(&mut m.triangles, inner_left, check_edge);
            check_vertex = apex(&m.triangles, check_edge);
        }

        while point(m.points, inner_right_apex)[1] < point(m.points, inner_right_org)[1] {
            next_self(inner_right);
            sym_self(&mut m.triangles, inner_right);
            inner_right_org = inner_right_apex;
            inner_right_apex = apex(&m.triangles, inner_right);
        }

        sym(&mut m.triangles, far_right, check_edge);
        check_vertex = apex(&m.triangles, check_edge);
        while point(m.points, check_vertex)[1] > point(m.points, far_right_pt)[1] {
            next(check_edge, far_right);
            far_right_pt = check_vertex;
            sym(&mut m.triangles, far_right, check_edge);
            check_vertex = apex(&m.triangles, check_edge);
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
    let left_cand = bump.alloc(HEdge::default());
    let right_cand = bump.alloc(HEdge::default());
    let base_edge = bump.alloc(HEdge::default());
    sym(&mut m.triangles, inner_left, left_cand);
    sym(&mut m.triangles, inner_right, right_cand);
    // Create the bottom new bounding triangle
    make_triangle(&mut m.triangles, base_edge);
    // Connect it to the bounding boxes of the left and right triangulations.
    bond(&mut m.triangles, base_edge, inner_left);
    next_self(base_edge);
    bond(&mut m.triangles, base_edge, inner_right);
    next_self(base_edge);
    set_org(&mut m.triangles, base_edge, inner_right_org);
    set_dest(&mut m.triangles, base_edge, inner_left_dest);

    // Fix the extreme triangles if necessary
    far_left_pt = org(&m.triangles, far_left);
    if inner_left_dest == far_left_pt {
        next(base_edge, far_left);
    }
    far_right_pt = dest(&m.triangles, far_right);
    if inner_right_org == far_right_pt {
        prev(base_edge, far_right);
    }

    // The vertices of the current knitting edge
    let mut lower_left = inner_left_dest;
    let mut lower_right = inner_right_org;
    // The candidate vertices for knitting
    let mut upper_left = apex(&m.triangles, left_cand);
    let mut upper_right = apex(&m.triangles, right_cand);
    // Walk up the gap between the two triangulations, knitting them together
    loop {
        let left_finished =
            counterclockwise(m.points, upper_left, lower_left, lower_right, bump) <= 0.0;
        let right_finished =
            counterclockwise(m.points, upper_right, lower_left, lower_right, bump) <= 0.0;
        let check_edge = bump.alloc(HEdge::default());
        let next_edge = bump.alloc(HEdge::default());
        if left_finished && right_finished {
            // Create the top new bounding triangle
            make_triangle(&mut m.triangles, next_edge);
            set_org(&mut m.triangles, next_edge, lower_left);
            set_dest(&mut m.triangles, next_edge, lower_right);
            // Apex is intentionally left INVALID
            // Connect it to the bounding boxes of the two triangulations
            bond(&mut m.triangles, next_edge, base_edge);
            next_self(next_edge);
            bond(&mut m.triangles, next_edge, right_cand);
            next_self(next_edge);
            bond(&mut m.triangles, next_edge, left_cand);

            // Special treatment for horizontal cuts
            if axis == 1 {
                far_left_pt = org(&m.triangles, far_left);
                far_right_pt = dest(&m.triangles, far_right);
                far_right_apex = apex(&m.triangles, far_right);
                sym(&mut m.triangles, far_left, check_edge);
                let mut check_vertex = apex(&m.triangles, check_edge);
                // The pointers to the extremal vertices are restored to the
                // leftmost and rightmost vertices (rather than topmost and
                // bottommost)
                while point(m.points, check_vertex)[0] < point(m.points, far_left_pt)[0] {
                    prev(check_edge, far_left);
                    far_left_pt = check_vertex;
                    sym(&mut m.triangles, far_left, check_edge);
                    check_vertex = apex(&m.triangles, check_edge);
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
            prev(left_cand, next_edge);
            sym_self(&mut m.triangles, next_edge);
            let mut next_apex = apex(&m.triangles, next_edge);
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
                let top_casing = bump.alloc(HEdge::default());
                let side_casing = bump.alloc(HEdge::default());
                let outer_casing = bump.alloc(HEdge::default());
                while bad_edge {
                    // Eliminate the edge with an edge flip.  As a result, the
                    // left triangulation will have one more boundary triangle.
                    next_self(next_edge);
                    sym(&mut m.triangles, next_edge, top_casing);
                    next_self(next_edge);
                    sym(&mut m.triangles, next_edge, side_casing);
                    bond(&mut m.triangles, next_edge, top_casing);
                    bond(&mut m.triangles, left_cand, side_casing);
                    next_self(left_cand);
                    sym(&mut m.triangles, left_cand, outer_casing);
                    prev_self(next_edge);
                    bond(&mut m.triangles, next_edge, outer_casing);
                    // Correct the vertices to reflect the edge flip
                    set_org(&mut m.triangles, left_cand, lower_left);
                    set_dest(&mut m.triangles, left_cand, INVALID_IND);
                    set_apex(&mut m.triangles, left_cand, next_apex);
                    set_org(&mut m.triangles, next_edge, INVALID_IND);
                    set_dest(&mut m.triangles, next_edge, upper_left);
                    set_apex(&mut m.triangles, next_edge, next_apex);
                    // Consider the newly exposed vertex
                    upper_left = next_apex;
                    // What vertex would be exposed if another edge were deleted?
                    copy(side_casing, next_edge);
                    next_apex = apex(&m.triangles, next_edge);
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
            next(right_cand, next_edge);
            sym_self(&mut m.triangles, next_edge);
            let mut next_apex = apex(&m.triangles, next_edge);
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
                let top_casing = bump.alloc(HEdge::default());
                let side_casing = bump.alloc(HEdge::default());
                let outer_casing = bump.alloc(HEdge::default());
                while bad_edge {
                    prev_self(next_edge);
                    sym(&mut m.triangles, next_edge, top_casing);
                    prev_self(next_edge);
                    sym(&mut m.triangles, next_edge, side_casing);
                    bond(&mut m.triangles, next_edge, top_casing);
                    bond(&mut m.triangles, right_cand, side_casing);
                    prev_self(right_cand);
                    sym(&mut m.triangles, right_cand, outer_casing);
                    next_self(next_edge);
                    bond(&mut m.triangles, next_edge, outer_casing);

                    set_org(&mut m.triangles, right_cand, INVALID_IND);
                    set_dest(&mut m.triangles, right_cand, lower_right);
                    set_apex(&mut m.triangles, right_cand, next_apex);
                    set_org(&mut m.triangles, next_edge, upper_right);
                    set_dest(&mut m.triangles, next_edge, INVALID_IND);
                    set_apex(&mut m.triangles, next_edge, next_apex);

                    upper_right = next_apex;

                    copy(side_casing, next_edge);
                    next_apex = apex(&m.triangles, next_edge);
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
            bond(&mut m.triangles, base_edge, right_cand);
            prev(right_cand, base_edge);
            set_dest(&mut m.triangles, base_edge, lower_left);
            lower_right = upper_right;
            sym(&mut m.triangles, base_edge, right_cand);
            upper_right = apex(&m.triangles, right_cand);
        } else {
            // Knit the triangulations, adding an edge from `upperleft'
            // to `lowerright'
            bond(&mut m.triangles, base_edge, left_cand);
            next(left_cand, base_edge);
            set_org(&mut m.triangles, base_edge, lower_right);
            lower_left = upper_left;
            sym(&mut m.triangles, base_edge, left_cand);
            upper_left = apex(&m.triangles, left_cand);
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
        let midtri = bump.alloc(HEdge::default());
        let tri1 = bump.alloc(HEdge::default());
        let tri2 = bump.alloc(HEdge::default());
        let tri3 = bump.alloc(HEdge::default());

        make_triangle(&mut m.triangles, midtri);
        make_triangle(&mut m.triangles, tri1);
        make_triangle(&mut m.triangles, tri2);
        make_triangle(&mut m.triangles, tri3);
        let area = counterclockwise(
            m.points,
            sorted_pt_inds[0],
            sorted_pt_inds[1],
            sorted_pt_inds[2],
            bump,
        );
        if area == 0.0 {
            // Three collinear vertices; the triangulation is two edges
            set_org(&mut m.triangles, midtri, sorted_pt_inds[0]);
            set_dest(&mut m.triangles, midtri, sorted_pt_inds[1]);
            set_org(&mut m.triangles, tri1, sorted_pt_inds[1]);
            set_dest(&mut m.triangles, tri1, sorted_pt_inds[0]);
            set_org(&mut m.triangles, tri2, sorted_pt_inds[2]);
            set_dest(&mut m.triangles, tri2, sorted_pt_inds[1]);
            set_org(&mut m.triangles, tri3, sorted_pt_inds[1]);
            set_dest(&mut m.triangles, tri3, sorted_pt_inds[2]);
            // All apices are intentionally left INVALID.
            bond(&mut m.triangles, midtri, tri1);
            bond(&mut m.triangles, tri2, tri3);
            next_self(midtri);
            prev_self(tri1);
            next_self(tri2);
            prev_self(tri3);
            bond(&mut m.triangles, midtri, tri3);
            bond(&mut m.triangles, tri1, tri2);
            next_self(midtri);
            prev_self(tri1);
            next_self(tri2);
            prev_self(tri3);
            bond(&mut m.triangles, midtri, tri1);
            bond(&mut m.triangles, tri2, tri3);

            copy(tri1, far_left);
            copy(tri2, far_right);
        } else {
            // The three vertices are not collinear; the triangulation is one
            // triangle, namely `midtri'.
            set_org(&mut m.triangles, midtri, sorted_pt_inds[0]);
            set_dest(&mut m.triangles, tri1, sorted_pt_inds[0]);
            set_org(&mut m.triangles, tri3, sorted_pt_inds[0]);
            if area > 0.0 {
                // The vertices are in counterclockwise order
                set_dest(&mut m.triangles, midtri, sorted_pt_inds[1]);
                set_org(&mut m.triangles, tri1, sorted_pt_inds[1]);
                set_dest(&mut m.triangles, tri2, sorted_pt_inds[1]);
                set_apex(&mut m.triangles, midtri, sorted_pt_inds[2]);
                set_org(&mut m.triangles, tri2, sorted_pt_inds[2]);
                set_dest(&mut m.triangles, tri3, sorted_pt_inds[2]);
            } else {
                // The vertices are in clockwise order
                set_dest(&mut m.triangles, midtri, sorted_pt_inds[2]);
                set_org(&mut m.triangles, tri1, sorted_pt_inds[2]);
                set_dest(&mut m.triangles, tri2, sorted_pt_inds[2]);
                set_apex(&mut m.triangles, midtri, sorted_pt_inds[1]);
                set_org(&mut m.triangles, tri2, sorted_pt_inds[1]);
                set_dest(&mut m.triangles, tri3, sorted_pt_inds[1]);
            }
            // The topology does not depend on how the vertices are ordered
            bond(&mut m.triangles, midtri, tri1);
            next_self(midtri);
            bond(&mut m.triangles, midtri, tri2);
            next_self(midtri);
            bond(&mut m.triangles, midtri, tri3);
            prev_self(tri1);
            next_self(tri2);
            bond(&mut m.triangles, tri1, tri2);
            prev_self(tri1);
            prev_self(tri3);
            bond(&mut m.triangles, tri1, tri3);
            next_self(tri2);
            prev_self(tri3);
            bond(&mut m.triangles, tri2, tri3);
            // Ensure that the origin of `farleft' is start
            copy(tri1, far_left);
            // Ensure that the destination of `farright' is `start + 2`
            if area > 0.0 {
                copy(tri2, far_right);
            } else {
                next(far_left, far_right);
            }
        }
    } else {
        let divider = len >> 1;
        let innerleft = bump.alloc(HEdge::default());
        let innerright = bump.alloc(HEdge::default());

        let (left, right) = sorted_pt_inds.split_at(divider);

        div_conq_recurse(m, left, 1 - axis, far_left, innerleft, bump);
        div_conq_recurse(m, right, 1 - axis, innerright, far_right, bump);

        merge_hulls(m, axis, far_left, innerleft, innerright, far_right, bump);
    }
}
