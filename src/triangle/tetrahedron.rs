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
    const n: usize = 3;
    const N: usize = 8;
    const mask: usize = 7;

    let mut i = 0;
    while i < N {
        gc[i] = i ^ (i >> 1);
        i += 1;
    }

    let mut e = 0;
    while e < N {
        let mut d = 0;
        while d < n {
            let f = e ^ (1 << d);
            let travel_bit = e ^ f;
            i = 0;
            while i < N {
                let k = gc[i] * (travel_bit << 1);
                let g = (k | (k / N)) & mask;
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
    const n: usize = 3;
    const N: usize = 8;
    let mut i = 1;
    while i < N {
        let mut v = !i; // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        let mut c = 0;
        while v != 0 {
            v >>= 1;
            c += 1;
        }
        tsb_1mod3[i] = c % n;
        i += 1;
    }
    tsb_1mod3
}

const enexttbl: [usize; 12] = enext_table();
const eprevtbl: [usize; 12] = eprev_table();
const esymtbl: [usize; 12] = [9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2];
const enextsymtbl: [usize; 12] = enext_esym_table(&esymtbl, &enexttbl);
const erepvesymtbl: [usize; 12] = eprev_esym_table(&esymtbl, &eprevtbl);
const bondtbl: [[usize; 12]; 12] = bond_table();
const fsymtbl: [[usize; 12]; 12] = fsym_table();
const facepivot1: [usize; 12] = face_pivot1(&esymtbl);
const facepivot2: [[usize; 12]; 12] = face_pivot2(&fsymtbl, &esymtbl);
const orgpivot: [usize; 12] = [3, 3, 1, 1, 2, 0, 0, 2, 1, 2, 3, 0];
const destpivot: [usize; 12] = [2, 0, 0, 2, 1, 2, 3, 0, 3, 3, 1, 1];
const apexpivot: [usize; 12] = [1, 2, 3, 0, 3, 3, 1, 1, 2, 0, 0, 2];
const oppopivot: [usize; 12] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
const vpivot: [usize; 4] = [11, 8, 9, 10];
const transgc: [[[usize; 8]; 3]; 8] = gray_code();
const tsb1mod3: [usize; 8] = trailing_set_bits_mod3();

#[inline(always)]
fn point(points: &[f64], idx: usize) -> &[f64] {
    &points[(idx * 3)..]
}

fn hilbert_split(
    points: &[f64],
    indices: &mut [usize],
    n_points: usize,
    gc0: usize,
    gc1: usize,
    bbox: &[f64],
) -> usize {
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
