#[derive(Clone)]
pub struct BBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl AsRef<BBox> for BBox {
    fn as_ref(&self) -> &BBox {
        self
    }
}

impl BBox {
    pub fn new(minx: f64, miny: f64, minz: f64, maxx: f64, maxy: f64, maxz: f64) -> BBox {
        BBox {
            min: [minx, miny, minz],
            max: [maxx, maxy, maxz],
        }
    }

    pub fn from_boxes<A: AsRef<BBox>, T: IntoIterator<Item = A>>(boxes: T) -> BBox {
        let mut res = BBox::default();
        for bbox in boxes {
            res.merge(bbox.as_ref());
        }
        res
    }

    #[inline]
    pub fn extend(&mut self, p: &[f64]) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(p[i]);
            self.max[i] = self.max[i].max(p[i]);
        }
    }

    pub fn merge(&mut self, other: &BBox) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(other.min[i]);
            self.max[i] = self.max[i].max(other.max[i]);
        }
    }

    pub fn scale(&mut self, s: f64) {
        let center = [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ];
        for i in 0..3 {
            self.min[i] = center[i] + s * (self.min[i] - center[i]);
            self.max[i] = center[i] + s * (self.max[i] - center[i]);
        }
    }

    pub fn scaled(mut self, s: f64) -> Self {
        self.scale(s);
        self
    }
}

impl Default for BBox {
    fn default() -> BBox {
        BBox {
            min: [f64::INFINITY, f64::INFINITY, f64::INFINITY],
            max: [f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY],
        }
    }
}

impl<A: AsRef<[f64]>> FromIterator<A> for BBox {
    fn from_iter<T: IntoIterator<Item = A>>(iter: T) -> Self {
        let mut bbox = BBox::default();
        for p in iter {
            bbox.extend(p.as_ref());
        }
        bbox
    }
}
