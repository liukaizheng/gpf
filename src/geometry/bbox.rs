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
    /// Check if this bounding box completely contains another bounding box.
    ///
    /// A box contains another if all points of the other box are within or on
    /// the boundary of this box. This includes the case where the boxes are identical.
    ///
    /// # Arguments
    /// * `other` - The bounding box to test for containment
    ///
    /// # Returns
    /// `true` if this box contains the other box, `false` otherwise
    ///
    /// # Examples
    /// ```
    /// use gpf::geometry::BBox;
    ///
    /// let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
    /// let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
    ///
    /// assert!(outer.contains(&inner));
    /// assert!(!inner.contains(&outer));
    /// ```
    pub fn contains(&self, other: &BBox) -> bool {
        for i in 0..3 {
            if self.min[i] > other.min[i] || self.max[i] < other.max[i] {
                return false;
            }
        }
        true
    }
    /// Check if this bounding box intersects with another bounding box.
    ///
    /// Two boxes intersect if they overlap in any way, including sharing faces,
    /// edges, or corners. Containment is also considered intersection.
    ///
    /// # Arguments
    /// * `other` - The bounding box to test for intersection
    ///
    /// # Returns
    /// `true` if the boxes intersect in any way, `false` if they are completely separate
    ///
    /// # Examples
    /// ```
    /// use gpf::geometry::BBox;
    ///
    /// let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
    /// let box_b = BBox::new(3.0, 3.0, 3.0, 8.0, 8.0, 8.0);  // Overlapping
    /// let box_c = BBox::new(10.0, 10.0, 10.0, 15.0, 15.0, 15.0);  // Separate
    ///
    /// assert!(box_a.intersects(&box_b));
    /// assert!(!box_a.intersects(&box_c));
    /// ```
    pub fn intersects(&self, other: &BBox) -> bool {
        for i in 0..3 {
            if self.max[i] < other.min[i] || self.min[i] > other.max[i] {
                return false;
            }
        }
        true
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
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_contains() {
        let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
        let partial = BBox::new(5.0, 5.0, 5.0, 15.0, 15.0, 15.0);

        assert!(outer.contains(&inner));
        assert!(!inner.contains(&outer));
        assert!(!outer.contains(&partial));
        assert!(!partial.contains(&outer));

        // Test identical boxes
        assert!(outer.contains(&outer));

        // Test edge case - box touching boundary
        let edge_box = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        assert!(outer.contains(&edge_box));
    }
    #[test]
    fn test_intersects() {
        let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
        let box_b = BBox::new(3.0, 3.0, 3.0, 8.0, 8.0, 8.0);
        let box_c = BBox::new(10.0, 10.0, 10.0, 15.0, 15.0, 15.0);

        // Overlapping boxes
        assert!(box_a.intersects(&box_b));
        assert!(box_b.intersects(&box_a));

        // Non-overlapping boxes
        assert!(!box_a.intersects(&box_c));
        assert!(!box_c.intersects(&box_a));

        // Box with itself
        assert!(box_a.intersects(&box_a));

        // Boxes touching at edge
        let touching = BBox::new(5.0, 5.0, 5.0, 10.0, 10.0, 10.0);
        assert!(box_a.intersects(&touching));

        // Boxes touching at corner
        let corner_touch = BBox::new(5.0, 5.0, 5.0, 8.0, 8.0, 8.0);
        assert!(box_a.intersects(&corner_touch));
    }
}
