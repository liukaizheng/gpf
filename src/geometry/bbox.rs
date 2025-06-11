#[derive(Clone)]
pub struct BBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

/// Represents the relationship between two bounding boxes.
/// 
/// This enum describes all possible spatial relationships between two 3D bounding boxes:
/// - Containment (one box completely inside another)
/// - Intersection (boxes overlap but neither contains the other)
/// - Disjoint (boxes don't touch or overlap at all)
/// 
/// # Examples
/// 
/// ```
/// use gpf::geometry::{BBox, BBoxRelation};
/// 
/// let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
/// let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
/// 
/// match outer.compare(&inner) {
///     BBoxRelation::AContainsB => println!("Outer box contains inner box"),
///     BBoxRelation::BContainsA => println!("Inner box contains outer box"),
///     BBoxRelation::Intersects => println!("Boxes intersect"),
///     BBoxRelation::Disjoint => println!("Boxes don't touch"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BBoxRelation {
    /// The first box (A) completely contains the second box (B).
    /// This includes the case where the boxes are identical.
    AContainsB,
    /// The second box (B) completely contains the first box (A).
    BContainsA,
    /// The boxes intersect or touch, but neither completely contains the other.
    /// This includes cases where boxes share faces, edges, or corners.
    Intersects,
    /// The boxes are completely separate with no intersection or contact.
    Disjoint,
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

    /// Compare this bounding box with another and return their spatial relationship.
    /// 
    /// This function determines one of four possible relationships between two bounding boxes:
    /// 1. This box contains the other box
    /// 2. The other box contains this box  
    /// 3. The boxes intersect but neither contains the other
    /// 4. The boxes are completely separate (disjoint)
    /// 
    /// # Arguments
    /// * `other` - The other bounding box to compare with
    /// 
    /// # Returns
    /// A `BBoxRelation` enum value indicating the spatial relationship:
    /// * `BBoxRelation::AContainsB` - if this box contains the other (including identical boxes)
    /// * `BBoxRelation::BContainsA` - if the other box contains this box
    /// * `BBoxRelation::Intersects` - if boxes intersect but neither contains the other
    /// * `BBoxRelation::Disjoint` - if boxes don't intersect at all
    /// 
    /// # Examples
    /// ```
    /// use gpf::geometry::{BBox, BBoxRelation};
    /// 
    /// // Container relationship
    /// let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
    /// let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
    /// assert_eq!(outer.compare(&inner), BBoxRelation::AContainsB);
    /// 
    /// // Intersection without containment
    /// let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
    /// let box_b = BBox::new(3.0, 3.0, 3.0, 8.0, 8.0, 8.0);
    /// assert_eq!(box_a.compare(&box_b), BBoxRelation::Intersects);
    /// 
    /// // Disjoint boxes
    /// let box_c = BBox::new(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    /// let box_d = BBox::new(5.0, 5.0, 5.0, 8.0, 8.0, 8.0);
    /// assert_eq!(box_c.compare(&box_d), BBoxRelation::Disjoint);
    /// ```
    pub fn compare(&self, other: &BBox) -> BBoxRelation {
        let self_contains_other = self.contains(other);
        let other_contains_self = other.contains(self);
        
        if self_contains_other && other_contains_self {
            // Boxes are identical or one is degenerate
            BBoxRelation::AContainsB
        } else if self_contains_other {
            BBoxRelation::AContainsB
        } else if other_contains_self {
            BBoxRelation::BContainsA
        } else if self.intersects(other) {
            BBoxRelation::Intersects
        } else {
            BBoxRelation::Disjoint
        }
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

    #[test]
    fn test_compare_a_contains_b() {
        let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
        
        assert_eq!(outer.compare(&inner), BBoxRelation::AContainsB);
    }

    #[test]
    fn test_compare_b_contains_a() {
        let inner = BBox::new(2.0, 2.0, 2.0, 8.0, 8.0, 8.0);
        let outer = BBox::new(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        
        assert_eq!(inner.compare(&outer), BBoxRelation::BContainsA);
    }

    #[test]
    fn test_compare_intersects() {
        let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
        let box_b = BBox::new(3.0, 3.0, 3.0, 8.0, 8.0, 8.0);
        
        assert_eq!(box_a.compare(&box_b), BBoxRelation::Intersects);
        assert_eq!(box_b.compare(&box_a), BBoxRelation::Intersects);
    }

    #[test]
    fn test_compare_disjoint() {
        let box_a = BBox::new(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let box_b = BBox::new(5.0, 5.0, 5.0, 8.0, 8.0, 8.0);
        
        assert_eq!(box_a.compare(&box_b), BBoxRelation::Disjoint);
        assert_eq!(box_b.compare(&box_a), BBoxRelation::Disjoint);
    }

    #[test]
    fn test_compare_identical() {
        let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
        let box_b = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
        
        assert_eq!(box_a.compare(&box_b), BBoxRelation::AContainsB);
        assert_eq!(box_b.compare(&box_a), BBoxRelation::AContainsB);
    }

    #[test]
    fn test_compare_edge_cases() {
        // Boxes touching at face
        let box_a = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 5.0);
        let box_b = BBox::new(5.0, 0.0, 0.0, 10.0, 5.0, 5.0);
        assert_eq!(box_a.compare(&box_b), BBoxRelation::Intersects);
        
        // Boxes touching at edge
        let box_c = BBox::new(5.0, 5.0, 0.0, 10.0, 10.0, 5.0);
        assert_eq!(box_a.compare(&box_c), BBoxRelation::Intersects);
        
        // Boxes touching at corner
        let box_d = BBox::new(5.0, 5.0, 5.0, 10.0, 10.0, 10.0);
        assert_eq!(box_a.compare(&box_d), BBoxRelation::Intersects);
        
        // Degenerate box (point)
        let point = BBox::new(2.5, 2.5, 2.5, 2.5, 2.5, 2.5);
        assert_eq!(box_a.compare(&point), BBoxRelation::AContainsB);
        
        // Degenerate box (line)
        let line = BBox::new(2.5, 2.5, 0.0, 2.5, 2.5, 5.0);
        assert_eq!(box_a.compare(&line), BBoxRelation::AContainsB);
        
        // Degenerate box (plane)
        let plane = BBox::new(0.0, 0.0, 2.5, 5.0, 5.0, 2.5);
        assert_eq!(box_a.compare(&plane), BBoxRelation::AContainsB);
    }

    #[test]
    fn test_compare_2d_cases() {
        // Test boxes that are flat in one dimension
        let flat_xy = BBox::new(0.0, 0.0, 0.0, 5.0, 5.0, 0.0);
        let flat_xy2 = BBox::new(-1.0, -1.0, 0.0, 6.0, 6.0, 0.0);
        assert_eq!(flat_xy.compare(&flat_xy2), BBoxRelation::BContainsA);
        
        let flat_intersect = BBox::new(3.0, 3.0, 0.0, 8.0, 8.0, 0.0);
        assert_eq!(flat_xy.compare(&flat_intersect), BBoxRelation::Intersects);
    }
}
