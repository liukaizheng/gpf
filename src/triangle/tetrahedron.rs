use std::collections::HashMap;

use crate::{
    math::{cross, norm, sub},
    predicates::{self, double_to_sign, sign_reverse, Orientation},
    INVALID_IND,
};
use bumpalo::{collections::Vec, vec, Bump};
use rand::seq::SliceRandom;

#[derive(Clone)]
pub struct TriFace {
    pub tet: usize,
    pub ver: usize, // version
}

impl Default for TriFace {
    #[inline(always)]
    fn default() -> Self {
        Self {
            tet: INVALID_IND,
            ver: 11,
        }
    }
}

impl TriFace {
    #[inline(always)]
    pub fn new(tet: usize, ver: usize) -> Self {
        Self { tet, ver }
    }

    #[inline(always)]
    pub fn set(&mut self, t: usize, v: usize) {
        self.tet = t;
        self.ver = v;
    }

    #[inline(always)]
    pub fn eprev_self(&mut self) {
        self.ver = EPREV_TBL[self.ver];
    }

    #[inline(always)]
    pub fn enext_self(&mut self) {
        self.ver = ENEXT_TBL[self.ver];
    }

    #[inline(always)]
    pub fn esym_self(&mut self) {
        self.ver = ESYM_TBL[self.ver];
    }

    #[inline(always)]
    pub fn eprev_esym_self(&mut self) {
        self.ver = EPREV_ESYM_TBL[self.ver];
    }

    #[inline(always)]
    pub fn enext_esym_self(&mut self) {
        self.ver = ENEXT_ESYM_TBL[self.ver];
    }
}

pub struct Tet {
    pub data: [usize; 4],
    pub nei: [TriFace; 4],
    mask: usize,
}

impl Default for Tet {
    #[inline(always)]
    fn default() -> Self {
        Self {
            data: [INVALID_IND; 4],
            nei: [
                TriFace::default(),
                TriFace::default(),
                TriFace::default(),
                TriFace::default(),
            ],
            mask: 0,
        }
    }
}

impl Tet {
    #[inline(always)]
    fn new(pa: usize, pb: usize, pc: usize, pd: usize) -> Self {
        Self {
            data: [pa, pb, pc, pd],
            nei: [
                TriFace::default(),
                TriFace::default(),
                TriFace::default(),
                TriFace::default(),
            ],
            mask: 0,
        }
    }

    #[inline(always)]
    pub fn index(&self, vid: usize) -> Option<usize> {
        self.data.iter().position(|&id| id == vid)
    }
}

pub struct TetMesh<'a, 'b: 'a> {
    pub points: &'a [f64],
    pub n_points: usize,
    pub tets: Vec<'b, Tet>,
    pub p2t: Vec<'b, usize>,
}

impl<'a, 'b: 'a> TetMesh<'a, 'b> {
    fn new(points: &'a [f64], bump: &'b Bump) -> Self {
        let n_points = points.len() / 3;
        Self {
            points,
            n_points,
            tets: Vec::new_in(bump),
            p2t: vec![in bump; INVALID_IND; n_points + 1],
        }
    }

    #[inline(always)]
    pub fn org(&self, f: &TriFace) -> usize {
        self.tets[f.tet].data[ORG_PIVOT[f.ver]]
    }

    #[inline(always)]
    pub fn dest(&self, f: &TriFace) -> usize {
        self.tets[f.tet].data[DEST_PIVOT[f.ver]]
    }

    #[inline(always)]
    pub fn apex(&self, f: &TriFace) -> usize {
        self.tets[f.tet].data[APEX_PIVOT[f.ver]]
    }

    #[inline(always)]
    pub fn oppo(&self, f: &TriFace) -> usize {
        self.tets[f.tet].data[OPPO_PIVOT[f.ver]]
    }

    #[inline(always)]
    pub fn infect(&mut self, t: usize) {
        self.tets[t].mask |= 1;
    }
    #[inline(always)]
    pub fn uninfect(&mut self, t: usize) {
        self.tets[t].mask &= !1;
    }
    #[inline(always)]
    pub fn infected(&self, t: usize) -> bool {
        (self.tets[t].mask & 1) != 0
    }

    #[inline(always)]
    pub fn mark_test(&mut self, t: usize) {
        self.tets[t].mask |= 2;
    }
    #[inline(always)]
    pub fn unmark_test(&mut self, t: usize) {
        self.tets[t].mask &= !2;
    }
    #[inline(always)]
    pub fn mark_tested(&self, t: usize) -> bool {
        (self.tets[t].mask & 2) != 0
    }

    #[inline(always)]
    fn bond(&mut self, t1: &TriFace, t2: &TriFace) {
        let f1 = &mut self.tets[t1.tet].nei[t1.ver & 3];
        f1.tet = t2.tet;
        f1.ver = BOND_TBL[t1.ver][t2.ver];

        let f2 = &mut self.tets[t2.tet].nei[t2.ver & 3];
        f2.tet = t1.tet;
        f2.ver = BOND_TBL[t2.ver][t1.ver];
    }

    #[inline(always)]
    pub fn fsym(&self, t1: &TriFace, t2: &mut TriFace) {
        let nf = &self.tets[t1.tet].nei[t1.ver & 3];
        t2.tet = nf.tet;
        t2.ver = FSYM_TBL[t1.ver][nf.ver];
    }

    #[inline(always)]
    pub fn fsym_self(&self, t: &mut TriFace) {
        let nf = &self.tets[t.tet].nei[t.ver & 3];
        t.tet = nf.tet;
        t.ver = FSYM_TBL[t.ver][nf.ver];
    }

    #[inline(always)]
    pub fn fnext_self(&self, t: &mut TriFace) {
        let nf = &self.tets[t.tet].nei[FACE_PIVOT1[t.ver]];
        t.tet = nf.tet;
        t.ver = FACE_PIVOT2[t.ver][nf.ver];
    }

    #[inline(always)]
    pub fn is_hull_tet(&self, idx: usize) -> bool {
        self.tets[idx].data[3] == self.n_points
    }

    #[inline(always)]
    pub fn point(&self, idx: usize) -> &[f64] {
        point(self.points, idx)
    }

    #[inline(always)]
    pub fn orient3d(&self, pa: usize, pb: usize, pc: usize, pd: usize) -> f64 {
        predicates::orient3d(
            &self.points[(pa * 3)..],
            &self.points[(pb * 3)..],
            &self.points[(pc * 3)..],
            &self.points[(pd * 3)..],
            self.tets.bump(),
        )
    }

    #[inline(always)]
    pub fn incident(&mut self, vid: usize) -> Vec<'b, usize> {
        let bump = self.tets.bump();
        let mut result = bumpalo::vec![in bump; self.p2t[vid]];
        self.mark_test(result[0]);
        let mut idx = 0;
        while idx < result.len() {
            let t = result[idx];
            let pos = self.tets[t].index(vid).unwrap();
            for i in 0..4 {
                if i == pos {
                    continue;
                }
                let nei = self.tets[t].nei[i].tet;
                if !self.mark_tested(nei) && !self.is_hull_tet(nei) {
                    self.mark_test(nei);
                    result.push(nei);
                }
            }
            idx += 1;
        }
        for &tid in &result {
            self.unmark_test(tid);
        }
        result
    }
}

const fn enext_table() -> [usize; 12] {
    let mut ret = [0usize; 12];
    let mut i = 0;
    while i < 12 {
        ret[i] = (i + 4) % 12;
        i += 1;
    }
    ret
}

const fn eprev_table() -> [usize; 12] {
    let mut ret = [0usize; 12];
    let mut i = 0;
    while i < 12 {
        ret[i] = (i + 8) % 12;
        i += 1;
    }
    ret
}

const fn enext_esym_table(esym_tbl: &[usize], enext_tbl: &[usize]) -> [usize; 12] {
    let mut ret = [0usize; 12];
    let mut i = 0;
    while i < 12 {
        ret[i] = esym_tbl[enext_tbl[i]];
        i += 1;
    }
    ret
}

const fn eprev_esym_table(esym_tbl: &[usize], eprev_tbl: &[usize]) -> [usize; 12] {
    let mut ret = [0usize; 12];
    let mut i = 0;
    while i < 12 {
        ret[i] = esym_tbl[eprev_tbl[i]];
        i += 1;
    }
    ret
}

const fn bond_table() -> [[usize; 12]; 12] {
    let mut ret = [[0usize; 12]; 12];
    let mut i = 0;
    while i < 12 {
        let mut j = 0;
        while j < 12 {
            ret[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
            j += 1;
        }
        i += 1;
    }
    ret
}

const fn fsym_table() -> [[usize; 12]; 12] {
    let mut ret = [[0usize; 12]; 12];
    let mut i = 0;
    while i < 12 {
        let mut j = 0;
        while j < 12 {
            ret[i][j] = (j + 12 - (i & 12)) % 12;
            j += 1;
        }
        i += 1;
    }
    ret
}

const fn face_pivot1(esym_tbl: &[usize]) -> [usize; 12] {
    let mut ret = [0usize; 12];
    let mut i = 0;
    while i < 12 {
        ret[i] = esym_tbl[i] & 3;
        i += 1;
    }
    ret
}

const fn face_pivot2(fsym_tbl: &[[usize; 12]], esym_tbl: &[usize]) -> [[usize; 12]; 12] {
    let mut ret = [[0usize; 12]; 12];
    let mut i = 0;
    while i < 12 {
        let mut j = 0;
        while j < 12 {
            ret[i][j] = fsym_tbl[esym_tbl[i]][j];
            j += 1;
        }
        i += 1;
    }
    ret
}

/// "Compact Hilbert Indices", Technical Report CS-2006-07
const fn gray_code() -> [[[usize; 8]; 3]; 8] {
    let mut trans_gc = [[[0usize; 8]; 3]; 8];
    let mut gc = [0usize; 8];
    const N1: usize = 3;
    const N2: usize = 8;
    const MASK: usize = 7;

    let mut i = 0;
    while i < N2 {
        gc[i] = i ^ (i >> 1);
        i += 1;
    }

    let mut e = 0;
    while e < N2 {
        let mut d = 0;
        while d < N1 {
            let f = e ^ (1 << d);
            let travel_bit = e ^ f;
            i = 0;
            while i < N2 {
                let k = gc[i] * (travel_bit << 1);
                let g = (k | (k / N2)) & MASK;
                trans_gc[e][d][i] = g ^ e;
                i += 1;
            }
            d += 1;
        }
        e += 1;
    }
    return trans_gc;
}

const fn trailing_set_bits_mod3() -> [usize; 8] {
    let mut tsb_1mod3 = [0usize; 8];
    tsb_1mod3[0] = 0;
    const N1: usize = 3;
    const N2: usize = 8;
    let mut i = 1;
    while i < N2 {
        let mut v = !i; // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        let mut c = 0;
        while v != 0 {
            v >>= 1;
            c += 1;
        }
        tsb_1mod3[i] = c % N1;
        i += 1;
    }
    tsb_1mod3
}

const ENEXT_TBL: [usize; 12] = enext_table();
const EPREV_TBL: [usize; 12] = eprev_table();
const ESYM_TBL: [usize; 12] = [9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2];
const ENEXT_ESYM_TBL: [usize; 12] = enext_esym_table(&ESYM_TBL, &ENEXT_TBL);
const EPREV_ESYM_TBL: [usize; 12] = eprev_esym_table(&ESYM_TBL, &EPREV_TBL);
const BOND_TBL: [[usize; 12]; 12] = bond_table();
const FSYM_TBL: [[usize; 12]; 12] = fsym_table();
const FACE_PIVOT1: [usize; 12] = face_pivot1(&ESYM_TBL);
const FACE_PIVOT2: [[usize; 12]; 12] = face_pivot2(&FSYM_TBL, &ESYM_TBL);
const ORG_PIVOT: [usize; 12] = [3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0];
const DEST_PIVOT: [usize; 12] = [2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1];
const APEX_PIVOT: [usize; 12] = [1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2];
const OPPO_PIVOT: [usize; 12] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
const TRANS_GC: [[[usize; 8]; 3]; 8] = gray_code();
const TSB1_MOD3: [usize; 8] = trailing_set_bits_mod3();
pub const VPIVOT: [usize; 4] = [11, 8, 9, 10];
const EPIVOT: [usize; 12] = [4, 5, 2, 11, 4, 5, 2, 11, 4, 5, 2, 11];

/// prev edge
#[inline(always)]
fn eprev(t1: &TriFace, t2: &mut TriFace) {
    t2.tet = t1.tet;
    t2.ver = EPREV_TBL[t1.ver];
}

/// next edge
fn enext(t1: &TriFace, t2: &mut TriFace) {
    t2.tet = t1.tet;
    t2.ver = ENEXT_TBL[t1.ver];
}

/// symmetric edge
#[inline(always)]
fn esym(t1: &TriFace, t2: &mut TriFace) {
    t2.tet = t1.tet;
    t2.ver = ESYM_TBL[t1.ver];
}

/// prev sym edge
#[inline(always)]
fn eprev_esym(t1: &TriFace, t2: &mut TriFace) {
    t2.tet = t1.tet;
    t2.ver = EPREV_ESYM_TBL[t1.ver];
}

/// next sym edge
#[inline(always)]
fn enext_esym(t1: &TriFace, t2: &mut TriFace) {
    t2.tet = t1.tet;
    t2.ver = ENEXT_ESYM_TBL[t1.ver];
}

#[inline(always)]
fn point(points: &[f64], idx: usize) -> &[f64] {
    let start = idx * 3;
    &points[start..(start + 3)]
}

fn hilbert_split(
    points: &[f64],
    indices: &mut [usize],
    gc0: usize,
    gc1: usize,
    bbox: &[f64],
) -> usize {
    let n_points = indices.len();
    if n_points == 1 {
        return 0;
    }
    let axis = (gc0 ^ gc1) >> 1;
    let split = (bbox[axis] + bbox[axis + 3]) * 0.5;

    let d = if (gc0 & (1 << axis)) == 0 { 1 } else { -1 };
    let mut i = 0;
    let mut j = 0;
    if d > 0 {
        loop {
            while i < n_points {
                if point(points, indices[i])[axis] >= split {
                    break;
                }
                i += 1;
            }
            while j < n_points {
                if point(points, indices[n_points - 1 - j])[axis] < split {
                    break;
                }
                j += 1;
            }
            if i + j == n_points {
                break;
            }
            indices.swap(i, n_points - 1 - j);
        }
    } else {
        loop {
            while i < n_points {
                if point(points, indices[i])[axis] <= split {
                    break;
                }
                i += 1;
            }
            while j < n_points {
                if point(points, indices[n_points - 1 - j])[axis] > split {
                    break;
                }
                j += 1;
            }
            if i + j == n_points {
                break;
            }
            indices.swap(i, n_points - 1 - j);
        }
    }
    i
}

struct SortOption {
    threshold: usize,
    hilbert_order: usize,
    hilbert_limit: usize,
    ratio: f64,
}

fn hilbert_sort(
    points: &[f64],
    indices: &mut [usize],
    bbox: &[f64],
    option: &SortOption,
    e: usize,
    d: usize,
    depth: usize,
) {
    let mid = hilbert_split(points, indices, TRANS_GC[e][d][3], TRANS_GC[e][d][4], bbox);
    let (left, right) = indices.split_at_mut(mid);

    let l_mid = hilbert_split(points, left, TRANS_GC[e][d][1], TRANS_GC[e][d][2], bbox);
    let (l_left, l_right) = left.split_at_mut(l_mid);
    let ll_mid = hilbert_split(points, l_left, TRANS_GC[e][d][0], TRANS_GC[e][d][1], bbox);
    let rl_mid = hilbert_split(points, l_right, TRANS_GC[e][d][2], TRANS_GC[e][d][3], bbox);

    let r_mid = hilbert_split(points, right, TRANS_GC[e][d][5], TRANS_GC[e][d][6], bbox);
    let (r_left, r_right) = right.split_at_mut(r_mid);
    let lr_mid = hilbert_split(points, r_left, TRANS_GC[e][d][4], TRANS_GC[e][d][5], bbox);
    let rr_mid = hilbert_split(points, r_right, TRANS_GC[e][d][6], TRANS_GC[e][d][7], bbox);

    if option.hilbert_order > 0 {
        if depth + 1 == option.hilbert_order {
            return;
        }
    }
    let (ll_left, ll_right) = l_left.split_at_mut(ll_mid);
    let (rl_left, rl_right) = l_right.split_at_mut(rl_mid);
    let (lr_left, lr_right) = r_left.split_at_mut(lr_mid);
    let (rr_left, rr_right) = r_right.split_at_mut(rr_mid);
    let arr = [
        ll_left, ll_right, rl_left, rl_right, lr_left, lr_right, rr_left, rr_right,
    ];
    const MASK: usize = 7;
    const N: usize = 3;
    for w in 0..8 {
        if arr[w].len() > option.hilbert_limit {
            let mut e_w;
            let mut k;
            if w == 0 {
                e_w = 0;
            } else {
                //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
                k = 2 * ((w - 1) / 2);
                e_w = k ^ (k >> 1); // = gc(k).
            }
            k = e_w;
            e_w = ((k << (d + 1)) & MASK) | ((k >> (N - d - 1)) & MASK);
            let ei = e ^ e_w;
            let d_w;
            if w == 0 {
                d_w = 0;
            } else {
                d_w = if (w % 2) == 0 {
                    TSB1_MOD3[w - 1]
                } else {
                    TSB1_MOD3[w]
                };
            }
            let di = (d + d_w + 1) % N;
            // Calculate the bounding box of the sub-box.
            let mut sbox = [0.0; 6];
            if (TRANS_GC[e][d][w] & 1) != 0 {
                // x-axis
                sbox[0] = (bbox[0] + bbox[3]) * 0.5;
                sbox[3] = bbox[3];
            } else {
                sbox[0] = bbox[0];
                sbox[3] = (bbox[0] + bbox[3]) * 0.5;
            }
            if (TRANS_GC[e][d][w] & 2) != 0 {
                // y-axis
                sbox[1] = (bbox[1] + bbox[4]) * 0.5;
                sbox[4] = bbox[4];
            } else {
                sbox[1] = bbox[1];
                sbox[4] = (bbox[1] + bbox[4]) * 0.5;
            }
            if (TRANS_GC[e][d][w] & 4) != 0 {
                // z-axis
                sbox[2] = (bbox[2] + bbox[5]) * 0.5;
                sbox[5] = bbox[5];
            } else {
                sbox[2] = bbox[2];
                sbox[5] = (bbox[2] + bbox[5]) * 0.5;
            }
            hilbert_sort(points, arr[w], &sbox, option, ei, di, depth + 1);
        }
    }
}

fn brio_multiscale_sort(
    points: &[f64],
    indices: &mut [usize],
    bbox: &[f64],
    option: &SortOption,
    depth: &mut usize,
) {
    let n_points = indices.len();
    if n_points >= option.threshold {
        *depth += 1;
        let mid = (n_points as f64 * option.ratio) as usize;
        let (left, right) = indices.split_at_mut(mid);
        brio_multiscale_sort(points, left, bbox, option, depth);
        hilbert_sort(points, right, bbox, option, 0, 0, 0);
    } else {
        hilbert_sort(points, indices, bbox, option, 0, 0, 0);
    }
}

#[inline(always)]
fn make_tet(tets: &mut TetMesh, face: &mut TriFace, pa: usize, pb: usize, pc: usize, pd: usize) {
    face.tet = tets.tets.len();
    face.ver = 11;
    tets.tets.push(Tet::new(pa, pb, pc, pd));
}

fn initial_delaunay(tets: &mut TetMesh, pa: usize, pb: usize, pc: usize, pd: usize) -> TriFace {
    let mut firsttet = TriFace::default();
    let mut tetopa = TriFace::default();
    let mut tetopb = TriFace::default();
    let mut tetopc = TriFace::default();
    let mut tetopd = TriFace::default();

    make_tet(tets, &mut firsttet, pa, pb, pc, pd);
    make_tet(tets, &mut tetopa, pb, pc, pd, tets.n_points);
    make_tet(tets, &mut tetopb, pc, pa, pd, tets.n_points);
    make_tet(tets, &mut tetopc, pa, pb, pd, tets.n_points);
    make_tet(tets, &mut tetopd, pb, pa, pc, tets.n_points);

    // Connect hull tetrahedra to firsttet (at four faces of firsttet).
    let mut worktet = TriFace::default();
    tets.bond(&firsttet, &tetopd);
    esym(&firsttet, &mut worktet);
    tets.bond(&worktet, &mut tetopc); // ab
    enext_esym(&firsttet, &mut worktet);
    tets.bond(&worktet, &tetopa); // bc
    eprev_esym(&firsttet, &mut worktet);
    tets.bond(&worktet, &tetopb); // ca

    // Connect hull tetrahedra together (at six edges of firsttet).
    let mut worktet1 = TriFace::default();
    esym(&tetopc, &mut worktet);
    esym(&tetopd, &mut worktet1);
    tets.bond(&worktet, &worktet1); // ab
    esym(&tetopa, &mut worktet);
    eprev_esym(&tetopd, &mut worktet1);
    tets.bond(&worktet, &worktet1); // bc
    esym(&tetopb, &mut worktet);
    enext_esym(&tetopd, &mut worktet1);
    tets.bond(&worktet, &worktet1); // ca
    eprev_esym(&tetopc, &mut worktet);
    enext_esym(&tetopb, &mut worktet1);
    tets.bond(&worktet, &worktet1); // da
    eprev_esym(&tetopa, &mut worktet);
    enext_esym(&tetopc, &mut worktet1);
    tets.bond(&worktet, &worktet1); // db
    eprev_esym(&tetopb, &mut worktet);
    enext_esym(&tetopa, &mut worktet1);
    tets.bond(&worktet, &worktet1); // dc

    tets.p2t[pa] = 0;
    tets.p2t[pb] = 0;
    tets.p2t[pc] = 0;
    tets.p2t[pd] = 0;

    tets.p2t[tets.n_points] = tetopa.tet;

    firsttet
}

enum LocateResult {
    OUTSIDE,
    ONVERTEX,
    ONEDGE,
    ONFACE,
    INTETRAHEDRON,
}

fn locate_dt(tets: &TetMesh, pid: usize, searchtet: &mut TriFace) -> LocateResult {
    if tets.is_hull_tet(searchtet.tet) {
        searchtet.tet = tets.tets[searchtet.tet].nei[3].tet;
    }

    searchtet.ver = 0;
    while searchtet.ver < 4 {
        let ori = tets.orient3d(
            tets.org(searchtet),
            tets.dest(searchtet),
            tets.apex(searchtet),
            pid,
        );
        if ori < 0.0 {
            break;
        }
        searchtet.ver += 1;
    }

    let loc;
    loop {
        let toppo = tets.oppo(searchtet);

        // Check if the vertex is we seek.
        if toppo == pid {
            // Adjust the origin of searchtet to be searchpt.
            searchtet.esym_self();
            searchtet.eprev_self();
            loc = LocateResult::ONVERTEX; // return ONVERTEX;
            break;
        }

        let oriorg = tets.orient3d(tets.dest(searchtet), tets.apex(searchtet), toppo, pid);
        if oriorg < 0.0 {
            searchtet.enext_esym_self();
        } else {
            let oridest = tets.orient3d(tets.apex(searchtet), tets.org(searchtet), toppo, pid);
            if oridest < 0.0 {
                searchtet.eprev_esym_self();
            } else {
                let oriapex = tets.orient3d(tets.org(searchtet), tets.dest(searchtet), toppo, pid);
                if oriapex < 0.0 {
                    searchtet.esym_self();
                } else {
                    // oriorg >= 0, oridest >= 0, oriapex >= 0 ==> found the point.
                    // The point we seek must be on the boundary of or inside this
                    //   tetrahedron. Check for boundary cases first.
                    if oriorg == 0.0 {
                        // Go to the face opposite to origin.
                        searchtet.enext_esym_self();
                        if oridest == 0.0 {
                            searchtet.eprev_self(); // edge oppo->apex
                            if oriapex == 0.0 {
                                // oppo is duplicated with p.
                                loc = LocateResult::ONVERTEX; // return ONVERTEX;
                                break;
                            }
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        if oriapex == 0.0 {
                            searchtet.enext_self(); // edge dest->oppo
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if oridest == 0.0 {
                        // Go to the face opposite to destination.
                        searchtet.eprev_esym_self();
                        if oriapex == 0.0 {
                            searchtet.eprev_self(); // edge oppo->org
                            loc = LocateResult::ONEDGE; // return ONEDGE;
                            break;
                        }
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    if oriapex == 0.0 {
                        // Go to the face opposite to apex
                        searchtet.esym_self();
                        loc = LocateResult::ONFACE; // return ONFACE;
                        break;
                    }
                    loc = LocateResult::INTETRAHEDRON;
                    break;
                }
            }
        }

        let nf = &tets.tets[searchtet.tet].nei[searchtet.ver & 3];
        searchtet.tet = nf.tet;
        searchtet.ver = nf.ver;

        if tets.is_hull_tet(searchtet.tet) {
            loc = LocateResult::OUTSIDE;
            break;
        }
    }
    loc
}

// Insphere test with symbolic perturbation
fn insphere_s(
    tets: &TetMesh,
    pa: usize,
    pb: usize,
    pc: usize,
    pd: usize,
    pe: usize,
) -> Orientation {
    let points = tets.points;
    let sign = predicates::insphere(
        &points[(pa * 3)..],
        &points[(pb * 3)..],
        &points[(pc * 3)..],
        &points[(pd * 3)..],
        &points[(pe * 3)..],
        tets.tets.bump(),
    );
    if sign != 0.0 {
        return double_to_sign(sign);
    }

    let mut pt = [pa, pb, pc, pd, pe];

    let mut swaps = 0; // Record the total number of swaps.
    let mut n = 5;
    while n > 0 {
        let mut count = 0;
        n = n - 1;
        let mut i = 0;
        while i < n {
            if pt[i] > pt[i + 1] {
                pt.swap(i, i + 1);
                count += 1;
            }
            i += 1;
        }
        swaps += count;
        // Continue if some points are swapped.
        if count == 0 {
            break;
        }
    }

    let mut ori = tets.orient3d(pt[1], pt[2], pt[3], pt[4]);
    if ori == 0.0 {
        ori = -tets.orient3d(pt[0], pt[2], pt[3], pt[4]);
    }

    if (swaps & 1) != 0 {
        sign_reverse(double_to_sign(ori))
    } else {
        double_to_sign(ori)
    }
}

// Insert a vertex using the Bowyer-Watson algorithm
fn insert_vertex_bw(
    tets: &mut TetMesh,
    pid: usize,
    searchtet: &mut TriFace,
    bw_face_map: &mut HashMap<(usize, usize), TriFace>,
) -> bool {
    let bump = tets.tets.bump();
    let mut cave_oldtet_list = Vec::new_in(bump);

    match locate_dt(tets, pid, searchtet) {
        LocateResult::OUTSIDE | LocateResult::INTETRAHEDRON => {
            tets.infect(searchtet.tet);
            cave_oldtet_list.push(searchtet.tet);
        }
        LocateResult::ONVERTEX => {
            return false;
        }
        LocateResult::ONEDGE => {
            let mut spintet = searchtet.clone();
            loop {
                tets.infect(spintet.tet);
                cave_oldtet_list.push(spintet.tet);
                tets.fnext_self(&mut spintet);
                if spintet.tet == searchtet.tet {
                    break;
                }
            }
        }
        LocateResult::ONFACE => {
            tets.infect(searchtet.tet);
            cave_oldtet_list.push(searchtet.tet);
            let nei_tet = tets.tets[searchtet.tet].nei[searchtet.ver & 3].tet;
            tets.infect(nei_tet);
            cave_oldtet_list.push(nei_tet);
        }
    }
    let mut cavetid;
    let mut neightid;
    let mut cave_bdry_list = Vec::new_in(bump);
    let mut idx = 0;
    while idx < cave_oldtet_list.len() {
        cavetid = cave_oldtet_list[idx];
        for ver in 0..4 {
            neightid = tets.tets[cavetid].nei[ver].tet;
            if tets.infected(neightid) {
                continue;
            }
            let mut enqflag = false;
            if !tets.mark_tested(neightid) {
                let pts = &tets.tets[neightid].data;
                if !tets.is_hull_tet(neightid) {
                    enqflag = insphere_s(tets, pts[0], pts[1], pts[2], pts[3], pid)
                        == Orientation::Negative;
                } else {
                    let ori = tets.orient3d(pts[0], pts[1], pts[2], pid);
                    if ori < 0.0 {
                        enqflag = true;
                    } else if ori == 0.0 {
                        let neineitet = tets.tets[neightid].nei[3].tet;
                        let nei_pts = &tets.tets[neineitet].data;
                        enqflag =
                            insphere_s(tets, nei_pts[0], nei_pts[1], nei_pts[2], nei_pts[3], pid)
                                == Orientation::Negative;
                    }
                }
                tets.mark_test(neightid);
            }

            if enqflag {
                tets.infect(neightid);
                cave_oldtet_list.push(neightid);
            } else {
                // A boundary face.
                cave_bdry_list.push(TriFace::new(cavetid, ver));
            }
        }
        idx += 1;
    }

    const ROW_V08_TBL: [usize; 12] = [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7];
    const ROW_V11_TBL: [usize; 12] = [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7];
    const COL_V01_TBL: [usize; 12] = [1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9];
    const COL_V02_TBL: [usize; 12] = [2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10];
    const COL_V08_TBL: [usize; 12] = [8, 8, 8, 8, 0, 0, 0, 0, 4, 4, 4, 4];
    const COL_V11_TBL: [usize; 12] = [11, 11, 11, 11, 3, 3, 3, 3, 7, 7, 7, 7];

    let f_out = cave_bdry_list.len();
    let v_out = (f_out + 4) >> 1;

    bw_face_map.clear();

    if v_out < 1024 {
        // pid to local vertex id
        for oldtet in &mut cave_bdry_list {
            // Get the tet outside the cavity.
            let mut neightet = tets.tets[oldtet.tet].nei[oldtet.ver].clone();
            tets.unmark_test(neightet.tet);

            if tets.is_hull_tet(oldtet.tet) {
                // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
                neightet.ver = EPIVOT[neightet.ver];
            }

            // Create a new tet in the cavity.
            let verts = [
                tets.dest(&neightet),
                tets.org(&neightet),
                tets.apex(&neightet),
            ];
            let mut newtet = TriFace::default();
            make_tet(tets, &mut newtet, verts[1], verts[0], pid, verts[2]);
            tets.tets[newtet.tet].nei[2] = neightet.clone();
            tets.tets[neightet.tet].nei[neightet.ver & 3]
                .set(newtet.tet, COL_V02_TBL[neightet.ver]);

            bw_face_map.insert((verts[0], verts[1]), TriFace::new(newtet.tet, 11));
            bw_face_map.insert((verts[1], verts[2]), TriFace::new(newtet.tet, 1));
            bw_face_map.insert((verts[2], verts[0]), TriFace::new(newtet.tet, 8));

            for vid in verts {
                let tid = tets.p2t[vid];
                if tets.tets[tid].data[2] != pid {
                    tets.p2t[vid] = newtet.tet;
                }
            }

            *oldtet = newtet;
        }

        // randomly pick a new tet
        *searchtet = cave_bdry_list[f_out >> 1].clone();
        tets.p2t[pid] = searchtet.tet;

        for mut neightet in cave_bdry_list {
            // TriFace neightet = cave_bdry_list[i];
            if tets.tets[neightet.tet].nei[3].tet == INVALID_IND {
                neightet.ver = 11;
                let j = tets.org(&neightet);
                let k = tets.dest(&neightet);
                let neineitet = bw_face_map.get(&(j, k)).unwrap();
                tets.tets[neightet.tet].nei[3].set(neineitet.tet, ROW_V11_TBL[neineitet.ver]);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3]
                    .set(neightet.tet, COL_V11_TBL[neineitet.ver]);
            }
            if tets.tets[neightet.tet].nei[1].tet == INVALID_IND {
                neightet.ver = 1;
                let j = tets.org(&neightet);
                let k = tets.dest(&neightet);
                let neineitet = bw_face_map.get(&(j, k)).unwrap();
                tets.tets[neightet.tet].nei[1].set(neineitet.tet, neineitet.ver);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3]
                    .set(neightet.tet, COL_V01_TBL[neineitet.ver]);
            }
            if tets.tets[neightet.tet].nei[0].tet == INVALID_IND {
                neightet.ver = 8;
                let j = tets.org(&neightet);
                let k = tets.dest(&neightet);
                let neineitet = bw_face_map.get(&(j, k)).unwrap();
                tets.tets[neightet.tet].nei[0].set(neineitet.tet, ROW_V08_TBL[neineitet.ver]);
                tets.tets[neineitet.tet].nei[neineitet.ver & 3]
                    .set(neightet.tet, COL_V08_TBL[neineitet.ver]);
            }
        }
    } else {
        // Fill a very large cavity with original neighboring searching method.
        for oldtet in &cave_bdry_list {
            let mut neightet = tets.tets[oldtet.tet].nei[oldtet.ver].clone();

            tets.unmark_test(neightet.tet);
            if tets.is_hull_tet(oldtet.tet) {
                // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
                neightet.ver = EPIVOT[neightet.ver];
            }

            let v = [
                tets.dest(&neightet),
                tets.org(&neightet),
                tets.apex(&neightet),
            ];
            let mut newtet = TriFace::default();
            make_tet(tets, &mut newtet, v[1], v[0], pid, v[2]);

            tets.tets[newtet.tet].nei[2].set(neightet.tet, neightet.ver);
            tets.tets[neightet.tet].nei[neightet.ver & 3]
                .set(newtet.tet, COL_V02_TBL[neightet.ver]);

            for vid in v {
                let tid = tets.p2t[vid];
                if tets.tets[tid].data[2] != pid {
                    tets.p2t[vid] = newtet.tet;
                }
            }
        }

        // randomly pick a new tet
        *searchtet = cave_bdry_list[f_out >> 1].clone();
        tets.p2t[pid] = searchtet.tet;

        for mut oldtet in cave_bdry_list {
            let mut neightet = TriFace::default();
            let mut newtet = TriFace::default();
            tets.fsym(&oldtet, &mut neightet);
            tets.fsym(&neightet, &mut newtet);
            // Oldtet and newtet must be at the same directed edge.
            // Connect the three other faces of this newtet.
            for _ in 0..3 {
                esym(&newtet, &mut neightet); // Go to the face.
                if tets.tets[neightet.tet].nei[neightet.ver & 3].tet == INVALID_IND {
                    // Find the adjacent face of this newtet.
                    let mut spintet = oldtet.clone();
                    loop {
                        tets.fnext_self(&mut spintet);
                        if !tets.infected(spintet.tet) {
                            break;
                        }
                    }
                    let mut neineitet = TriFace::default();
                    tets.fsym(&spintet, &mut neineitet);
                    neineitet.esym_self();
                    tets.bond(&neightet, &neineitet);
                }
                newtet.enext_self();
                oldtet.enext_self();
            }
        }
    }

    for tid in cave_oldtet_list {
        tets.tets[tid].mask = !0;
    }

    return true;
}

pub fn tetrahedralize<'a, 'b: 'a>(points: &'a [f64], bump: &'b Bump) -> TetMesh<'a, 'b> {
    let cmp = |x: &&f64, y: &&f64| x.partial_cmp(y).unwrap();
    let bbox = [
        // min_corner
        *points.iter().step_by(3).min_by(cmp).unwrap(),
        *points[1..].iter().step_by(3).min_by(cmp).unwrap(),
        *points[2..].iter().step_by(3).min_by(cmp).unwrap(),
        // max_corner
        *points.iter().step_by(3).max_by(cmp).unwrap(),
        *points[1..].iter().step_by(3).max_by(cmp).unwrap(),
        *points[2..].iter().step_by(3).max_by(cmp).unwrap(),
    ];
    let mut mesh = TetMesh::new(points, bump);
    let mut sorted_pt_inds = Vec::from_iter_in(0..mesh.n_points, bump);
    sorted_pt_inds.shuffle(&mut rand::thread_rng());
    const SORT_OPTION: SortOption = SortOption {
        threshold: 64,
        hilbert_order: 52,
        hilbert_limit: 8,
        ratio: 0.125,
    };
    const EPS: f64 = 1e-6;
    let mut n_group = 1;
    brio_multiscale_sort(
        points,
        &mut sorted_pt_inds,
        &bbox,
        &SORT_OPTION,
        &mut n_group,
    );

    let mut temp_vec3 = vec![in bump; 0.0, 0.0, 0.0];
    let bbox_size1 = {
        sub(&bbox[3..6], &bbox, &mut temp_vec3);
        norm(&temp_vec3)
    };
    let bbox_size2 = bbox_size1 * bbox_size1;
    let bbox_size3 = bbox_size2 * bbox_size1;
    let mut i = 1;
    {
        let p0 = point(points, sorted_pt_inds[0]);
        while i < mesh.n_points {
            let pi = point(points, sorted_pt_inds[i]);
            sub(p0, pi, &mut temp_vec3);
            if norm(&temp_vec3) / bbox_size1 < EPS {
                i += 1;
            } else {
                break;
            }
        }
        if i > 1 {
            if i == mesh.n_points {
                return mesh;
            }
            sorted_pt_inds.swap(1, i);
        }

        i = 2;
        let p1 = point(points, sorted_pt_inds[1]);
        sub(p1, p0, &mut temp_vec3);
        while i < mesh.n_points {
            let pi = point(points, sorted_pt_inds[i]);
            let mut v2 = vec![in bump; 0.0; 3];
            sub(pi, p0, &mut v2);
            let mut n = vec![in bump; 0.0; 3];
            cross(&temp_vec3, &v2, &mut n);
            if norm(&n) / bbox_size2 < EPS {
                i += 1;
            } else {
                break;
            }
        }
        if i > 2 {
            if i == mesh.n_points {
                return mesh;
            };
            sorted_pt_inds.swap(2, i);
        }
    }
    i = 3;
    let mut ori = 0.0;
    while i < mesh.n_points {
        ori = predicates::orient3d_fast(
            point(points, sorted_pt_inds[0]),
            point(points, sorted_pt_inds[1]),
            point(points, sorted_pt_inds[2]),
            point(points, sorted_pt_inds[i]),
        );
        if ori.abs() / bbox_size3 < EPS {
            i += 1;
        } else {
            break;
        }
    }
    if i > 3 {
        if i == mesh.n_points {
            return mesh;
        }
        sorted_pt_inds.swap(3, i);
    }

    // follow the right rule
    if ori > 0.0 {
        sorted_pt_inds.swap(0, 1);
    }

    let mut search_tet = initial_delaunay(
        &mut mesh,
        sorted_pt_inds[0],
        sorted_pt_inds[1],
        sorted_pt_inds[2],
        sorted_pt_inds[3],
    );
    let mut bw_face_map = HashMap::new();
    for i in 4..mesh.n_points {
        if !insert_vertex_bw(
            &mut mesh,
            sorted_pt_inds[i],
            &mut search_tet,
            &mut bw_face_map,
        ) {
            break;
        }
    }

    let mut count = 0;
    let mut tet_map = vec![in bump; 0; mesh.tets.len()];
    for i in 0..tet_map.len() {
        let t = &mesh.tets[i];
        if t.mask != !0 {
            tet_map[i] = count;
            count += 1;
        }
    }

    mesh.tets.retain(|t| t.mask != !0);

    // update neighbors
    for tet in &mut mesh.tets {
        for nei in &mut tet.nei {
            nei.tet = tet_map[nei.tet];
        }
    }
    // make p2t[i] not be a ghost
    for i in 0..mesh.n_points {
        mesh.p2t[i] = tet_map[mesh.p2t[i]];
        if mesh.is_hull_tet(mesh.p2t[i]) {
            mesh.p2t[i] = mesh.tets[mesh.p2t[i]].nei[3].tet;
        }
    }

    mesh
}
