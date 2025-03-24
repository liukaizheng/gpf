use hashbrown::HashMap;
use itertools::Itertools;
use std::alloc::Allocator;
use std::ops::{Add, Index, Mul};
use std::slice::SliceIndex;

use bumpalo::Bump;

use crate::math::{dot, sub_in};
use crate::mesh::{ElementId, FaceId, HalfedgeId, ManifoldMesh, Mesh, VertexId};
use crate::predicates::{ImplicitPointSSI, incircle};
use crate::predicates::{Orientation, Point2D, orient2d_2d};
use crate::{INVALID_IND, is_negative, is_positive, point, predicates, twin_index};

struct Triangulation<'a, A: Allocator + Copy> {
    points: &'a [f64],
    alloc: A,
    mesh: ManifoldMesh<A>,
    sorted_vertices: Vec<VertexId, A>,
}
impl<'a, A: Allocator + Copy> Triangulation<'a, A> {
    fn triangulate(mut self, is_horizontal: bool) -> (HalfedgeId, ManifoldMesh<A>) {
        let n_points = self.points.len() >> 1;
        self.mesh.new_vertices(n_points);
        self.mesh.reserve_edges(n_points * 6);
        let [bdy_hid, _] = self.div_conq_recurse(0, n_points, is_horizontal);
        (bdy_hid, self.mesh)
    }

    fn div_conq_recurse(
        &mut self,
        start: usize,
        end: usize,
        is_horizontal: bool,
    ) -> [HalfedgeId; 2] {
        let len = end - start;
        match len {
            2 => {
                let va = self.sorted_vertices[start];
                let vb = self.sorted_vertices[start + 1];
                let [ha, hb, hc] = self.new_triangle(va, vb, VertexId::default());
                let twin_ha = self.mesh.he_twin(ha);
                let twin_hb = self.mesh.he_twin(hb);
                let twin_hc = self.mesh.he_twin(hc);
                self.mesh
                    .new_face_by_halfedges(&[twin_hc, twin_hb, twin_ha]);
                if self.cmp(va, vb, !is_horizontal).is_lt() {
                    [twin_hc, twin_hb]
                } else {
                    [hb, hc]
                }
            }
            3 => {
                let va = self.sorted_vertices[start];
                let vb = self.sorted_vertices[start + 1];
                let vc = self.sorted_vertices[start + 2];
                let area = self.counterclockwise(va, vb, vc);
                let order = area.partial_cmp(&0.0).unwrap();
                match order {
                    std::cmp::Ordering::Less | std::cmp::Ordering::Greater => {
                        let halfedges = if order == std::cmp::Ordering::Greater {
                            self.new_triangle(va, vb, vc)
                        } else {
                            self.new_triangle(va, vc, vb)
                        };
                        let side_halfedges = halfedges.map(|hid| {
                            let vid = self.mesh.he_from(hid);
                            self.mesh.new_edge_by_veritces(vid, VertexId::default())
                        });

                        for ((&ha, &hb_twin), hc_twin) in side_halfedges
                            .iter()
                            .circular_tuple_windows()
                            .zip(halfedges)
                        {
                            self.mesh.new_face_by_halfedges(&[
                                ha,
                                self.mesh.he_twin(hb_twin),
                                self.mesh.he_twin(hc_twin),
                            ]);
                        }

                        let (&min_v_hid, &max_v_hid) = side_halfedges
                            .iter()
                            .minmax_by(|&&h1, &&h2| {
                                self.cmp(
                                    self.mesh.he_from(h1),
                                    self.mesh.he_from(h2),
                                    !is_horizontal,
                                )
                            })
                            .into_option()
                            .unwrap();
                        [min_v_hid, self.mesh.he_twin(max_v_hid)]
                    }
                    std::cmp::Ordering::Equal => {
                        let ha = self.mesh.new_edge_by_veritces(va, vb);
                        let hb = self.mesh.new_edge_by_veritces(vb, vc);
                        let halfedges = [ha, hb, self.mesh.he_twin(hb), self.mesh.he_twin(ha)];
                        let side_halfedges = halfedges.map(|hid| {
                            let vid = self.mesh.he_from(hid);
                            self.mesh.new_edge_by_veritces(vid, VertexId::default())
                        });

                        for ((&ha, &hb_twin), hc_twin) in side_halfedges
                            .iter()
                            .circular_tuple_windows()
                            .zip(halfedges)
                        {
                            self.mesh.new_face_by_halfedges(&[
                                ha,
                                self.mesh.he_twin(hb_twin),
                                self.mesh.he_twin(hc_twin),
                            ]);
                        }

                        if self.cmp(va, vc, !is_horizontal).is_lt() {
                            [side_halfedges[0], self.mesh.he_twin(side_halfedges[2])]
                        } else {
                            [side_halfedges[2], self.mesh.he_twin(side_halfedges[0])]
                        }
                    }
                }
            }
            _ => {
                let mid = start + (len >> 1);
                let [_, lr_hid] = self.div_conq_recurse(start, mid, !is_horizontal);
                let [rl_hid, _] = self.div_conq_recurse(mid, end, !is_horizontal);
                self.merge_hulls(lr_hid, rl_hid, is_horizontal)
            }
        }
    }

    fn merge_hulls(
        &mut self,
        mut lr_hid: HalfedgeId,
        mut rl_hid: HalfedgeId,
        is_horizontal: bool,
    ) -> [HalfedgeId; 2] {
        let [mut lt_vid, mut lb_vid] = self.he_dest_apex(lr_hid); // left top, left bottom
        let [mut rb_vid, mut rt_vid] = self.he_apex_org(rl_hid); // right bottom, right top

        loop {
            let mut stop = true;
            if self.counterclockwise(lt_vid, lb_vid, rt_vid) > 0.0 {
                stop = false;
                lr_hid = self.mesh.he_prev_twin(lr_hid);
                lt_vid = lb_vid;
                lb_vid = self.mesh.he_to_to(lr_hid);
            }

            if self.counterclockwise(rb_vid, rt_vid, lt_vid) > 0.0 {
                stop = false;
                rl_hid = self.mesh.he_next_twin(rl_hid);
                rt_vid = rb_vid;
                rb_vid = self.mesh.he_to_to(rl_hid);
            }
            if stop {
                break;
            }
        }

        let mut left_hid = self.mesh.he_twin(lr_hid);
        let mut right_hid = self.mesh.he_twin(rl_hid);
        lb_vid = lt_vid;
        rb_vid = rt_vid;

        lt_vid = self.mesh.he_to_to(left_hid);
        rt_vid = self.mesh.he_to_to(right_hid);

        let mut bottom_hid = self.mesh.new_edge_by_veritces(rb_vid, lb_vid);
        {
            let ha = self.mesh.new_edge_by_veritces(lb_vid, VertexId::default());
            let hb = self.mesh.new_edge_by_veritces(VertexId::default(), rb_vid);
            self.mesh.new_face_by_halfedges(&[bottom_hid, ha, hb]);

            self.mesh.he_replace(lr_hid, self.mesh.he_twin(ha));
            self.mesh.he_replace(rl_hid, self.mesh.he_twin(hb));
        }
        let mut top_hid = self.mesh.he_twin(bottom_hid);

        loop {
            let left_finished = self.counterclockwise(lt_vid, lb_vid, rb_vid) <= 0.0;
            let right_finished = self.counterclockwise(rb_vid, rt_vid, lb_vid) <= 0.0;
            if left_finished && right_finished {
                break;
            }

            if !left_finished {
                let mut curr_hid = self.mesh.he_prev_twin(left_hid);
                loop {
                    let apex_vid = self.mesh.he_to_to(curr_hid);
                    if !apex_vid.valid() {
                        break;
                    }

                    if self.incircle(lb_vid, rb_vid, lt_vid, apex_vid) <= 0.0 {
                        break;
                    }

                    self.mesh.flip(curr_hid);
                    lt_vid = apex_vid;
                    curr_hid = self.mesh.he_next_twin(curr_hid);
                    debug_assert!(self.mesh.he_from(curr_hid) == lb_vid);
                    debug_assert!(self.mesh.he_to(curr_hid) == apex_vid);
                }
                left_hid = self.mesh.he_next(self.mesh.he_twin(curr_hid));
            }

            if !right_finished {
                let mut curr_hid = self.mesh.he_next_twin(right_hid);
                loop {
                    let apex_vid = self.mesh.he_to_to(curr_hid);
                    if !apex_vid.valid() {
                        break;
                    }

                    if self.incircle(lb_vid, rb_vid, rt_vid, apex_vid) <= 0.0 {
                        break;
                    }

                    self.mesh.flip(curr_hid);
                    rt_vid = apex_vid;
                    curr_hid = self.mesh.he_prev_twin(self.mesh.he_twin(curr_hid));
                    debug_assert!(self.mesh.he_from(curr_hid) == rt_vid);
                    debug_assert!(self.mesh.he_to(curr_hid) == rb_vid);
                }
                right_hid = self.mesh.he_prev(self.mesh.he_twin(curr_hid));
            }

            if left_finished
                || (!right_finished && self.incircle(lt_vid, lb_vid, rb_vid, rt_vid) > 0.0)
            {
                let prev_right_hid = self.mesh.he_prev(right_hid);
                self.mesh.he_replace(right_hid, top_hid);
                let new_hid = self.mesh.new_edge_by_veritces(rt_vid, lb_vid);
                self.mesh.he_replace(prev_right_hid, new_hid);

                top_hid = self.mesh.he_twin(new_hid);
                right_hid = self.mesh.he_twin(prev_right_hid);

                rb_vid = rt_vid;
                rt_vid = self.mesh.he_to_to(right_hid);
            } else {
                let next_left_hid = self.mesh.he_next(left_hid);
                self.mesh.he_replace(left_hid, top_hid);
                let new_hid = self.mesh.new_edge_by_veritces(rb_vid, lt_vid);
                self.mesh.he_replace(next_left_hid, new_hid);

                top_hid = self.mesh.he_twin(new_hid);
                left_hid = self.mesh.he_twin(next_left_hid);

                lb_vid = lt_vid;
                lt_vid = self.mesh.he_to_to(left_hid);
            }
        }

        self.mesh.new_face_by_halfedges(&[
            top_hid,
            self.mesh.he_twin(right_hid),
            self.mesh.he_twin(left_hid),
        ]);

        if is_horizontal {
            top_hid = self.rotate_prev(top_hid, |pa, pb| pa[1] < pb[1]);
            top_hid = self.rotate_next(top_hid, |pa, pb| pa[1] >= pb[1]);

            bottom_hid = self.rotate_next(bottom_hid, |pa, pb| pa[1] <= pb[1]);
            bottom_hid = self.rotate_prev(bottom_hid, |pa, pb| pa[1] > pb[1]);
            [self.mesh.he_next(bottom_hid), self.mesh.he_prev(top_hid)]
        } else {
            top_hid = self.rotate_next(top_hid, |pa, pb| pa[0] <= pb[0]);
            top_hid = self.rotate_prev(top_hid, |pa, pb| pa[0] > pb[0]);

            bottom_hid = self.rotate_prev(bottom_hid, |pa, pb| pa[0] < pb[0]);
            bottom_hid = self.rotate_next(bottom_hid, |pa, pb| pa[0] >= pb[0]);
            [self.mesh.he_next(top_hid), self.mesh.he_prev(bottom_hid)]
        }
    }

    fn new_triangle(&mut self, va: VertexId, vb: VertexId, vc: VertexId) -> [HalfedgeId; 3] {
        let ha = self.mesh.new_edge_by_veritces(va, vb);
        let hb = self.mesh.new_edge_by_veritces(vb, vc);
        let hc = self.mesh.new_edge_by_veritces(vc, va);
        let halfedges = [ha, hb, hc];
        self.mesh.new_face_by_halfedges(&halfedges);
        halfedges
    }

    #[inline]
    fn he_apex_org(&self, hid: HalfedgeId) -> [VertexId; 2] {
        let next_hid = self.mesh.he_next(hid);
        let nnext_hid = self.mesh.he_next(next_hid);
        [self.mesh.he_to(next_hid), self.mesh.he_to(nnext_hid)]
    }

    #[inline]
    fn he_dest_apex(&self, hid: HalfedgeId) -> [VertexId; 2] {
        let next_hid = self.mesh.he_next(hid);
        [self.mesh.he_to(hid), self.mesh.he_to(next_hid)]
    }

    fn cmp(&self, va: VertexId, vb: VertexId, is_horizontal: bool) -> std::cmp::Ordering {
        let pa = point::<2>(self.points, va.0);
        let pb = point::<2>(self.points, vb.0);
        if is_horizontal {
            pa.partial_cmp(&pb).unwrap()
        } else {
            [pa[1], -pa[0]].partial_cmp(&[pb[1], -pa[0]]).unwrap()
        }
    }

    fn rotate_prev(
        &self,
        mut hid: HalfedgeId,
        stop_fn: impl Fn(&[f64], &[f64]) -> bool,
    ) -> HalfedgeId {
        let [va, vb] = self.mesh.he_vertices(hid);
        let mut pa = point::<2>(self.points, va.0);
        let mut pb = point::<2>(self.points, vb.0);
        loop {
            if stop_fn(pa, pb) {
                break;
            }
            hid = self.mesh.he_prev(self.mesh.he_prev_twin(hid));
            pb = pa;
            pa = point::<2>(self.points, self.mesh.he_from(hid).0);
        }
        hid
    }

    fn rotate_next(
        &self,
        mut hid: HalfedgeId,
        stop_fn: impl Fn(&[f64], &[f64]) -> bool,
    ) -> HalfedgeId {
        let [va, vb] = self.mesh.he_vertices(hid);
        let mut pa = point::<2>(self.points, va.0);
        let mut pb = point::<2>(self.points, vb.0);
        loop {
            if stop_fn(pa, pb) {
                break;
            }
            hid = self.mesh.he_next(self.mesh.he_next_twin(hid));
            pa = pb;
            pb = point::<2>(self.points, self.mesh.he_to(hid).0);
        }
        hid
    }

    #[inline]
    fn counterclockwise(&self, va: VertexId, vb: VertexId, vc: VertexId) -> f64 {
        predicates::orient2d(
            point::<2>(self.points, va.0),
            point::<2>(self.points, vb.0),
            point::<2>(self.points, vc.0),
            self.alloc,
        )
    }

    #[inline]
    fn incircle(&self, va: VertexId, vb: VertexId, vc: VertexId, vd: VertexId) -> f64 {
        predicates::incircle(
            point::<2>(self.points, va.0),
            point::<2>(self.points, vb.0),
            point::<2>(self.points, vc.0),
            point::<2>(self.points, vd.0),
            self.alloc,
        )
    }
}

/// A constrained Delaunay triangulation.
struct CDT<'a, A: Allocator + Copy> {
    points: &'a [f64],
    segments: &'a [usize],
    point_indices: Vec<Point2D, A>,
    mesh: ManifoldMesh<A>,
    halfedge_marks: Vec<usize, A>,
    alloc: A,
}

impl<'a, A: Allocator + Copy> CDT<'a, A> {
    fn perform(&mut self) {
        for (idx, seg) in self.segments.chunks(2).enumerate() {
            let va = VertexId(seg[0]);
            let vb = VertexId(seg[1]);
            if va != vb {
                self.insert_segment(va, vb, idx << 1);
            }
        }
    }

    fn insert_segment(&mut self, mut va: VertexId, mut vb: VertexId, mark: usize) {
        let start_hid = self.scout_segment(self.mesh.v_halfedge(va), vb, mark);
        va = self.mesh.he_from(start_hid);
        if va == vb {
            return;
        }
        let end_hid = self.scout_segment(self.mesh.v_halfedge(vb), va, twin_index(mark));
        vb = self.mesh.he_from(end_hid);
        if vb == va {
            return;
        }

        self.constrain(start_hid, vb, mark);
    }

    fn scout_segment(&mut self, mut hid: HalfedgeId, vb: VertexId, mark: usize) -> HalfedgeId {
        let mut right_vid = self.mesh.he_to(hid);
        if right_vid == vb {
            self.set_edge_mark(hid, mark);
            return self.mesh.v_halfedge(vb);
        }

        let va = self.mesh.he_from(hid);
        let mut left_vid = self.mesh.he_to_to(hid);

        let mut left_ori = self.orient(vb, va, left_vid);
        let mut right_ori = self.orient(vb, va, right_vid);
        loop {
            if left_ori.is_pos() || right_ori.is_pos() {
                break;
            }
            hid = self.mesh.he_twin_next(hid);
            left_ori = right_ori;
            right_vid = self.mesh.he_to(hid);
            right_ori = self.orient(vb, va, right_vid);
        }

        loop {
            if !right_ori.is_neg() && !left_ori.is_pos() {
                break;
            }
            hid = self.mesh.he_prev_twin(hid);
            right_ori = left_ori;
            left_vid = self.mesh.he_to_to(hid);
            left_ori = self.orient(vb, va, left_vid);
        }

        let mut collinear = false;
        if left_ori.is_zero() {
            hid = self.mesh.he_prev_twin(hid);
            collinear = true;
        } else if right_ori.is_zero() {
            collinear = true;
        }

        if collinear {
            self.set_edge_mark(hid, mark);
            let end_vid = self.mesh.he_to(hid);
            if end_vid == vb {
                return self.mesh.v_halfedge(vb);
            } else {
                return self.scout_segment(self.mesh.v_halfedge(end_vid), vb, mark);
            }
        } else {
            let split_hid = self.mesh.he_next(hid);
            if self.halfedge_marks[split_hid] == INVALID_IND {
                return hid;
            } else {
                let [_, new_hid] = self.split_edge(split_hid, mark);
                return self.scout_segment(new_hid, vb, mark);
            }
        }
    }

    fn constrain(&mut self, bottom_right_hid: HalfedgeId, vb: VertexId, mark: usize) {
        let va = self.mesh.he_from(bottom_right_hid);
        let mut flip_hid = self.mesh.he_next(bottom_right_hid);
        self.mesh.flip(flip_hid);
        loop {
            let top_vid = self.mesh.he_from(flip_hid);
            if top_vid == vb {
                let fixup_hid = self.mesh.he_twin_next(flip_hid);
                self.fixup(flip_hid, false);
                self.fixup(fixup_hid, true);
                self.set_edge_mark(flip_hid, twin_index(mark));
                break;
            }

            let ori = self.orient(va, vb, top_vid);
            if ori.is_zero() {
                let fixup_hid = self.mesh.he_twin_next(flip_hid);
                self.fixup(flip_hid, false);
                self.fixup(fixup_hid, true);
                self.set_edge_mark(flip_hid, twin_index(mark));

                let hid = self.scout_segment(self.mesh.he_prev_twin(flip_hid), vb, mark);
                if self.mesh.he_from(hid) != vb {
                    self.constrain(hid, vb, mark);
                }
                break;
            } else {
                if ori.is_pos() {
                    let fixup_hid = self.mesh.he_twin_next(flip_hid);
                    self.fixup(fixup_hid, true);
                    flip_hid = self.mesh.he_prev(flip_hid);
                } else {
                    self.fixup(flip_hid, false);
                    flip_hid = self.mesh.he_twin_next(flip_hid);
                }

                if self.halfedge_marks[flip_hid] != INVALID_IND {
                    let [new_hid1, new_hid2] = self.split_edge(flip_hid, mark);
                    self.fixup(self.mesh.he_twin(new_hid1), false);
                    self.fixup(self.mesh.he_next(new_hid1), true);
                    let hid = self.scout_segment(new_hid2, vb, mark);
                    if self.mesh.he_from(hid) != vb {
                        self.constrain(hid, vb, mark);
                    }
                    break;
                } else {
                    self.mesh.flip(flip_hid);
                }
            }
        }
    }

    fn fixup(&mut self, hid: HalfedgeId, left_side: bool) {
        let flip_hid = self.mesh.he_next(hid);
        if self.halfedge_marks[flip_hid] != INVALID_IND {
            return;
        }
        let twin_flip_hid = self.mesh.he_twin(flip_hid);

        let bottom_vid = self.mesh.he_to_to(twin_flip_hid);
        if !bottom_vid.valid() {
            return;
        }

        let top_vid = self.mesh.he_from(hid);
        let [left_vid, right_vid] = self.mesh.he_vertices(flip_hid);

        if left_side {
            if !self.orient(top_vid, left_vid, bottom_vid).is_pos() {
                return;
            }
        } else {
            if !self.orient(bottom_vid, right_vid, top_vid).is_pos() {
                return;
            }
        }

        if self.orient(left_vid, bottom_vid, right_vid).is_pos() {
            if !self
                .incircle(left_vid, bottom_vid, right_vid, top_vid)
                .is_pos()
            {
                return;
            }
        }

        self.mesh.flip(flip_hid);
        self.fixup(hid, left_side);
        self.fixup(twin_flip_hid, left_side);
    }

    fn split_edge(&mut self, hid: HalfedgeId, input_mark: usize) -> [HalfedgeId; 2] {
        let twin_hid = self.mesh.he_twin(hid);
        let fid = self.mesh.he_face(hid);
        let twin_fid = self.mesh.he_face(twin_hid);
        let bottom_vid = self.mesh.he_to_to(hid);
        let top_vid = self.mesh.he_to_to(twin_hid);

        let left_vid = self.mesh.he_to(hid);

        let mark = self.halfedge_marks[hid];
        debug_assert!(mark != INVALID_IND);

        let new_vid = self.mesh.split_edge(self.mesh.he_edge(hid));

        let vab = &self.segments[((input_mark >> 1) << 1)..];
        let vcd = &self.segments[((mark >> 1) << 1)..];

        self.point_indices.push(Point2D::I(ImplicitPointSSI::new(
            vab[0], vab[1], vcd[0], vcd[1],
        )));
        let new_hid = self.mesh.v_halfedge(new_vid);
        self.new_edge_call_back(new_hid);
        if self.mesh.he_to(new_hid) == left_vid {
            self.set_edge_mark(new_hid, mark);
        } else {
            self.set_edge_mark(new_hid, twin_index(mark));
        }

        let new_hid1 = self.mesh.split_face(fid, bottom_vid, new_vid);
        self.new_edge_call_back(new_hid1);
        self.set_edge_mark(new_hid1, input_mark);

        let new_hid2 = self.mesh.split_face(twin_fid, new_vid, top_vid);
        self.new_edge_call_back(new_hid2);
        [new_hid1, new_hid2]
    }

    fn new_edge_call_back(&mut self, new_hid: HalfedgeId) {
        if new_hid.0 >= self.halfedge_marks.len() {
            self.halfedge_marks
                .resize(self.halfedge_marks.len() + 2, INVALID_IND);
        }
    }

    #[inline]
    fn set_edge_mark(&mut self, hid: HalfedgeId, mark: usize) {
        self.halfedge_marks[hid] = mark;
        self.halfedge_marks[self.mesh.he_twin(hid)] = twin_index(mark);
    }

    #[inline]
    fn orient(&self, va: VertexId, vb: VertexId, vc: VertexId) -> Orientation {
        orient2d_2d::orient2d(
            &self.point_indices[va],
            &self.point_indices[vb],
            &self.point_indices[vc],
            self.points,
            self.alloc,
        )
    }

    #[inline]
    fn incircle(&self, va: VertexId, vb: VertexId, vc: VertexId, vd: VertexId) -> Orientation {
        incircle::incircle(
            &self.point_indices[va],
            &self.point_indices[vb],
            &self.point_indices[vc],
            &self.point_indices[vd],
            self.points,
            self.alloc,
        )
    }

    #[inline]
    fn face_is_ghost(&self, fid: FaceId) -> bool {
        self.mesh.face(fid).vertices().any(|v| !v.valid())
    }

    fn extract_invalid_faces(&self) -> Vec<FaceId, A> {
        let mut result = Vec::with_capacity_in(self.mesh.n_faces_capacity(), self.alloc);
        let mut visited = Vec::with_capacity_in(self.mesh.n_faces_capacity(), self.alloc);
        visited.resize(self.mesh.n_faces_capacity(), false);
        for face in self.mesh.faces() {
            let fid = *face;
            if visited[fid] {
                continue;
            }
            visited[fid] = true;

            let mut queue = Vec::new_in(self.alloc);
            queue.push(fid);
            let mut keep = !self.face_is_ghost(fid);
            let mut idx = 0;
            while idx < queue.len() {
                let fid = queue[idx];
                idx += 1;
                for he in self.mesh.face(fid).halfedges() {
                    let hid = *he;
                    let mark = self.halfedge_marks[hid];
                    if mark == INVALID_IND {
                        let adj_fid = *he.twin().face();
                        if visited[adj_fid] {
                            continue;
                        }
                        visited[adj_fid] = true;
                        if !keep && self.face_is_ghost(adj_fid) {
                            keep = false;
                        }
                        queue.push(adj_fid);
                    } else {
                        if is_negative(self.halfedge_marks[hid]) {
                            keep = false;
                        }
                    }
                }
            }
            if keep {
                result.extend(queue);
            }
        }
        result
    }
}

fn get_triangulated_mesh<A: Allocator + Copy>(
    points: &[f64],
    is_horizontal: bool,
    alloc: A,
) -> (HalfedgeId, ManifoldMesh<A>) {
    let n_points = points.len() >> 1;
    let mut sorted_vertices = Vec::<VertexId, _>::with_capacity_in(n_points, alloc);
    sorted_vertices.extend((0..n_points).map(|idx| VertexId(idx)));
    sorted_vertices.sort_unstable_by(|&i, &j| {
        point::<2>(points, i.0)
            .partial_cmp(point::<2>(points, j.0))
            .unwrap()
    });
    alternate_axes(points, &mut sorted_vertices, is_horizontal);
    let triangulation = Triangulation {
        points,
        alloc,
        mesh: ManifoldMesh::new(([] as [[usize; 0]; 0]).into_iter(), alloc),
        sorted_vertices,
    };
    triangulation.triangulate(is_horizontal)
}

pub fn triangulate_points<A: Allocator + Copy>(
    points: &[f64],
    is_horizontal: bool,
    alloc: A,
) -> Vec<usize, A> {
    let (_, mesh) = get_triangulated_mesh(points, is_horizontal, alloc);
    let mut result = Vec::with_capacity_in(mesh.n_faces(), alloc);
    result.extend(
        mesh.faces()
            .filter_map(|f| {
                let fid = *f;
                let ha = mesh.f_halfedge(fid);
                let hb = mesh.he_next(ha);
                let hc = mesh.he_next(hb);
                let va = mesh.he_to(ha);
                let vb = mesh.he_to(hb);
                let vc = mesh.he_to(hc);
                if va.valid() && vb.valid() && vc.valid() {
                    Some([va.0, vb.0, vc.0])
                } else {
                    None
                }
            })
            .flatten(),
    );
    result
}

pub fn triangulate1<A: Allocator + Copy>(
    points: &[f64],
    segments: &[usize],
    is_horizontal: bool,
    alloc: A,
) -> Vec<usize, A> {
    let (bdy_hid, mut mesh) = get_triangulated_mesh(points, is_horizontal, alloc);
    set_boundary_vertex_halfedges(&mut mesh, bdy_hid);
    let mut halfedge_marks = Vec::<usize, A>::with_capacity_in(mesh.n_halfedges_capacity(), alloc);
    halfedge_marks.extend(std::iter::repeat(INVALID_IND).take(mesh.n_halfedges_capacity()));
    let mut point_indices = Vec::with_capacity_in(mesh.n_vertices_capacity(), alloc);
    point_indices.extend((0..mesh.n_vertices_capacity()).map(|idx| Point2D::E(idx)));
    let mut cdt = CDT {
        points,
        segments,
        point_indices,
        mesh,
        halfedge_marks,
        alloc,
    };
    cdt.perform();
    let valid_faces = cdt.extract_invalid_faces();

    let mut result = Vec::<usize, A>::with_capacity_in(valid_faces.len() * 3, alloc);
    for fid in valid_faces.into_iter() {
        for he in cdt.mesh.face(fid).halfedges() {
            result.push(**he.to());
        }
    }
    result
}

pub fn triangulate_with_new_points<A: Allocator + Copy>(
    points: &[f64],
    segments: &[usize],
    is_horizontal: bool,
    alloc: A,
) -> (Vec<f64, A>, Vec<usize, A>) {
    let n_old_points = points.len() >> 1;
    let (bdy_hid, mut mesh) = get_triangulated_mesh(points, is_horizontal, alloc);
    set_boundary_vertex_halfedges(&mut mesh, bdy_hid);
    let mut halfedge_marks = Vec::<usize, A>::with_capacity_in(mesh.n_halfedges_capacity(), alloc);
    halfedge_marks.extend(std::iter::repeat(INVALID_IND).take(mesh.n_halfedges_capacity()));
    let mut point_indices = Vec::with_capacity_in(mesh.n_vertices_capacity(), alloc);
    point_indices.extend((0..mesh.n_vertices_capacity()).map(|idx| Point2D::E(idx)));
    let mut cdt = CDT {
        points,
        segments,
        point_indices,
        mesh,
        halfedge_marks,
        alloc,
    };
    cdt.perform();
    let mut new_points =
        Vec::with_capacity_in((cdt.mesh.n_vertices_capacity() - n_old_points) << 1, alloc);
    for i in n_old_points..cdt.mesh.n_vertices_capacity() {
        new_points.extend_from_slice(&cdt.point_indices[i].to_explicit(&cdt.points));
    }

    let valid_faces = cdt.extract_invalid_faces();

    let mut triangles = Vec::<usize, A>::with_capacity_in(valid_faces.len() * 3, alloc);
    for fid in valid_faces.into_iter() {
        for he in cdt.mesh.face(fid).halfedges() {
            triangles.push(**he.to());
        }
    }
    (new_points, triangles)
}

fn set_boundary_vertex_halfedges<A: Allocator + Copy>(
    mesh: &mut ManifoldMesh<A>,
    first_hid: HalfedgeId,
) {
    let mut curr_hid = first_hid;
    loop {
        let prev_hid = mesh.he_prev(curr_hid);
        let va = mesh.he_to(prev_hid);
        mesh.set_v_halfedge(va, mesh.he_twin(prev_hid));
        curr_hid = mesh.he_next_twin(curr_hid);
        if curr_hid == first_hid {
            break;
        }
    }
}

pub fn triangulate_face_into_mesh<
    U: AsRef<[f64]>,
    T: IntoIterator<Item = U>,
    A: Allocator + Copy,
>(
    face_points: T,
    alloc: A,
) {
    let mut points = Vec::new_in(alloc);
    let mut segments = Vec::new_in(alloc);
    let mut start = 0;
    for polygon in face_points {
        let loop_points = polygon.as_ref();
        points.extend_from_slice(loop_points);
        let end = (loop_points.len() >> 1) + start;
        for (i, j) in (start..end).circular_tuple_windows() {
            segments.push(i);
            segments.push(j);
        }
        start = end;
    }
}

#[inline(always)]
fn point3(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

#[inline]
pub fn triangulate_polygon<A: Allocator + Copy>(
    points: &[f64],
    segments: &[usize],
    o: &[f64],
    x: &[f64],
    y: &[f64],
    bump: A,
) -> Vec<usize, A> {
    let [new_segments, new_to_ori_map] = unique_indices(segments, bump);
    let mut points_2d = Vec::new_in(bump);
    points_2d.extend(
        new_to_ori_map
            .iter()
            .map(|&idx| {
                let p = point3(points, idx);
                let mut v = std::vec::from_elem_in(0.0, 3, bump);
                sub_in(p, o, &mut v);
                [dot(&v, x), dot(&v, y)]
            })
            .flatten(),
    );
    let mut result = Vec::new_in(bump);
    result.extend(
        triangulate1(&points_2d, &new_segments, true, bump)
            .into_iter()
            .map(|idx| new_to_ori_map[idx]),
    );
    result
}

#[inline]
pub fn triangulate_polygon_soup(
    points: &[f64],
    edges: &[Vec<usize>],
    axes: &[f64],
) -> (Vec<usize>, Vec<usize>) {
    let mut triangles = Vec::new();
    let mut parents = Vec::new();
    let mut bump = Bump::new();
    for (idx, (segments, axis_data)) in edges.iter().zip(axes.chunks(9)).enumerate() {
        bump.reset();
        let face_triangles = triangulate_polygon(
            points,
            segments,
            &axis_data[0..3],
            &axis_data[3..6],
            &axis_data[6..9],
            &bump,
        );
        parents.resize(parents.len() + face_triangles.len() / 3, idx);
        triangles.extend(face_triangles);
    }
    (triangles, parents)
}

fn unique_indices<A: Allocator + Copy>(indices: &[usize], bump: A) -> [Vec<usize, A>; 2] {
    let mut count = 0;
    let mut map = HashMap::with_capacity_in(indices.len(), bump);
    let mut result = Vec::with_capacity_in(indices.len(), bump);
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
    let mut new_to_ori_map = std::vec::from_elem_in(0, map.len(), bump);
    for (k, v) in map {
        new_to_ori_map[v] = k;
    }
    [result, new_to_ori_map]
}

pub(crate) fn alternate_axes<T: Copy + Add<usize, Output = T> + Mul<usize, Output = T>>(
    points: &[f64],
    indices: &mut [T],
    mut is_horizontal: bool,
) where
    [f64]: Index<T, Output = f64>,
{
    let len = indices.len();
    let divider = len >> 1;
    if len <= 3 {
        is_horizontal = true;
    }

    if is_horizontal {
        indices.select_nth_unstable_by(divider, |&i, &j| {
            let i = i * 2;
            let j = j * 2;
            (points[i], points[i + 1])
                .partial_cmp(&(points[j], points[j + 1]))
                .unwrap()
        });
    } else {
        indices.select_nth_unstable_by(divider, |&i, &j| {
            let i = i * 2;
            let j = j * 2;
            (points[i + 1], -points[i])
                .partial_cmp(&(points[j + 1], -points[j]))
                .unwrap()
        });
    }

    let (left, right) = indices.split_at_mut(divider);
    if len - divider >= 2 {
        if divider >= 2 {
            alternate_axes(points, left, !is_horizontal);
        }
        alternate_axes(points, right, !is_horizontal);
    }
}
