use std::{
    alloc::Allocator,
    collections::VecDeque,
    ops::{Bound, RangeBounds},
};

use crate::INVALID_IND;

#[derive(PartialEq, Eq)]
enum Color {
    Red,
    Black,
    #[allow(dead_code)]
    Gray,
}

struct Node<V> {
    value: V,
    left: usize,
    right: usize,
    parent: usize,
    color: Color,
}

impl<V> Node<V> {
    fn new(value: V) -> Self {
        Node {
            value,
            left: INVALID_IND,
            right: INVALID_IND,
            parent: INVALID_IND,
            color: Color::Red,
        }
    }
}

pub struct RBTree<V, A: Allocator + Copy> {
    root: usize,
    nodes: Vec<Node<V>, A>,
    removed_nodes: Option<VecDeque<usize, A>>,
}

impl<V: PartialOrd, A: Allocator + Copy> RBTree<V, A> {
    pub fn into_iter(self) -> impl Iterator<Item = V> {
        self.nodes.into_iter().filter_map(|node| match node.color {
            Color::Gray => None,
            _ => Some(node.value),
        })
    }
    /// Returns an iterator over the values in the tree
    pub fn inorder_iter<'a>(&'a self) -> impl Iterator<Item = &'a V> {
        let mut stack = Vec::new();

        // Initialize stack with leftmost path from root
        if self.root != INVALID_IND {
            self.push_left_path(&mut stack, self.root);
        }

        std::iter::from_fn(move || {
            // Pop the next node from stack
            let node_id = stack.pop()?;
            let value = &self.nodes[node_id].value;

            // If this node has a right child, push all left nodes from right subtree
            if self.nodes[node_id].right != INVALID_IND {
                self.push_left_path(&mut stack, self.nodes[node_id].right);
            }

            Some(value)
        })
    }

    /// Returns an iterator over the values in the tree that are in the specified range
    pub fn range<'a, U, R>(&'a self, range: R) -> impl Iterator<Item = &'a V>
    where
        R: RangeBounds<U>,
        V: PartialOrd<U>,
        U: PartialOrd<V>,
    {
        let mut stack = Vec::new();

        // Build initial stack for range start
        if self.root != INVALID_IND {
            self.build_range_stack(&mut stack, self.root, &range);
        }

        std::iter::from_fn(move || {
            // Pop the next node from stack
            let node_id = stack.pop()?;
            let value = &self.nodes[node_id].value;

            // Check if we've reached the end bound
            match range.end_bound() {
                Bound::Included(bound_value) => {
                    if value.partial_cmp(bound_value).unwrap().is_gt() {
                        return None;
                    }
                }
                Bound::Excluded(bound_value) => {
                    if value.partial_cmp(bound_value).unwrap().is_ge() {
                        return None;
                    }
                }
                Bound::Unbounded => {}
            }

            // If this node has a right child, push all left nodes from right subtree
            if self.nodes[node_id].right != INVALID_IND {
                self.push_left_path(&mut stack, self.nodes[node_id].right);
            }

            Some(value)
        })
    }

    /// Returns a mutable iterator over the values in the tree that are in the specified range
    pub fn range_mut<'a, U, R>(&'a mut self, range: R) -> impl Iterator<Item = &'a mut V>
    where
        R: RangeBounds<U>,
        V: PartialOrd<U>,
        U: PartialOrd<V>,
    {
        let mut stack = Vec::new();

        // Build initial stack for range start
        if self.root != INVALID_IND {
            self.build_range_stack(&mut stack, self.root, &range);
        }

        // Get a raw pointer to the nodes for safe access within the closure
        let nodes_ptr = self.nodes.as_mut_ptr();
        let nodes_len = self.nodes.len();

        std::iter::from_fn(move || {
            // Pop the next node from stack
            let node_id = stack.pop()?;

            // Safety check: ensure node_id is within bounds
            if node_id >= nodes_len {
                return None;
            }

            // SAFETY: We know node_id is valid and within bounds
            let node_value = unsafe { &(*nodes_ptr.add(node_id)).value };

            // Check if we've reached the end bound before returning the value
            match range.end_bound() {
                Bound::Included(bound_value) => {
                    if node_value.partial_cmp(bound_value).unwrap().is_gt() {
                        return None;
                    }
                }
                Bound::Excluded(bound_value) => {
                    if node_value.partial_cmp(bound_value).unwrap().is_ge() {
                        return None;
                    }
                }
                Bound::Unbounded => {}
            }

            // If this node has a right child, push all left nodes from right subtree
            // SAFETY: We know node_id is valid and within bounds
            let right_child = unsafe { (*nodes_ptr.add(node_id)).right };
            if right_child != INVALID_IND {
                let mut current = right_child;
                while current != INVALID_IND && current < nodes_len {
                    stack.push(current);
                    // SAFETY: We know current is valid and within bounds
                    current = unsafe { (*nodes_ptr.add(current)).left };
                }
            }

            // SAFETY: We know node_id is valid and we're returning a mutable reference
            // with the same lifetime as the iterator, which ensures exclusive access
            unsafe { Some(&mut (*nodes_ptr.add(node_id)).value) }
        })
    }

    /// Helper function to push all nodes along the leftmost path starting from the given node
    fn push_left_path(&self, stack: &mut Vec<usize>, mut node_id: usize) {
        while node_id != INVALID_IND {
            stack.push(node_id);
            node_id = self.nodes[node_id].left;
        }
    }

    /// Helper function to build initial stack for range queries
    fn build_range_stack<U, R>(&self, stack: &mut Vec<usize>, mut current: usize, range: &R)
    where
        R: RangeBounds<U>,
        V: PartialOrd<U>,
        U: PartialOrd<V>,
    {
        // Navigate to the appropriate starting point and build stack
        while current != INVALID_IND {
            let curr_val = &self.nodes[current].value;

            let go_right = match range.start_bound() {
                Bound::Included(v) => curr_val.partial_cmp(v).unwrap().is_lt(),
                Bound::Excluded(v) => curr_val.partial_cmp(v).unwrap().is_le(),
                Bound::Unbounded => false,
            };

            if go_right {
                current = self.nodes[current].right;
            } else {
                stack.push(current);
                current = self.nodes[current].left;
            }
        }
    }
    pub fn new(alloc: A, reuse: bool) -> Self {
        RBTree {
            root: INVALID_IND,
            nodes: Vec::new_in(alloc),
            removed_nodes: if reuse {
                Some(VecDeque::new_in(alloc))
            } else {
                None
            },
        }
    }

    fn new_node(&mut self, value: V) -> usize {
        if let Some(removed_nodes) = &mut self.removed_nodes {
            if let Some(index) = removed_nodes.pop_back() {
                let node = &mut self.nodes[index];
                node.value = value;
                node.color = Color::Red;
                return index;
            }
        }
        let index = self.nodes.len();
        self.nodes.push(Node::new(value));
        index
    }

    pub fn insert(&mut self, val: V) {
        let new_node_id = self.new_node(val);
        self.bst_insert(new_node_id);
        self.insert_case1(new_node_id);
        println!("the root is {}", self.root);
        for (i, node) in self.nodes.iter().enumerate() {
            println!(
                "node: {}, parent: {}, left: {}, right: {}",
                i, node.parent, node.left, node.right
            );
        }
    }

    pub fn bst_insert(&mut self, new_node_id: usize) {
        let mut parent_id = INVALID_IND;
        {
            let v = &self.nodes[new_node_id].value;
            let mut curr_id = self.root;
            loop {
                if curr_id == INVALID_IND {
                    break;
                }
                parent_id = curr_id;
                let curr_node = &self.nodes[curr_id];
                if v.partial_cmp(&curr_node.value).unwrap().is_lt() {
                    curr_id = curr_node.left;
                } else {
                    curr_id = curr_node.right;
                }
            }
        }
        self.nodes[new_node_id].parent = parent_id;
        if parent_id != INVALID_IND {
            if self.nodes[new_node_id]
                .value
                .partial_cmp(&self.nodes[parent_id].value)
                .unwrap()
                .is_lt()
            {
                self.nodes[parent_id].left = new_node_id;
            } else {
                self.nodes[parent_id].right = new_node_id;
            }
        }
    }

    fn insert_case1(&mut self, node_id: usize) {
        if self.nodes[node_id].parent == INVALID_IND {
            self.nodes[node_id].color = Color::Black;
            self.root = node_id;
            if self.root == 1 {
                println!("bug here");
            }
        } else {
            self.insert_case2(node_id);
        }
    }

    fn insert_case2(&mut self, node_id: usize) {
        match self.nodes[self.nodes[node_id].parent].color {
            Color::Black => {}
            _ => {
                self.insert_case3(node_id);
            }
        }
    }

    fn insert_case3(&mut self, node_id: usize) {
        let parent_id = self.nodes[node_id].parent;
        let grandparent_id = self.nodes[parent_id].parent;
        let uncle_id = self.twin(parent_id, grandparent_id);
        if uncle_id != INVALID_IND && self.nodes[uncle_id].color == Color::Red {
            self.nodes[parent_id].color = Color::Black;
            self.nodes[uncle_id].color = Color::Black;
            self.nodes[grandparent_id].color = Color::Red;
            self.insert_case1(grandparent_id);
        } else {
            self.insert_case4(node_id);
        }
    }

    fn insert_case4(&mut self, mut node_id: usize) {
        let parent_id = self.nodes[node_id].parent;
        let grandparent_id = self.nodes[parent_id].parent;
        if node_id == self.nodes[parent_id].right && parent_id == self.nodes[grandparent_id].left {
            self.rotate_left(node_id);
            node_id = self.nodes[node_id].left;
        } else if node_id == self.nodes[parent_id].left
            && parent_id == self.nodes[grandparent_id].right
        {
            self.rotate_right(node_id);
            node_id = self.nodes[node_id].right;
        }
        self.insert_case5(node_id);
    }

    fn insert_case5(&mut self, node_id: usize) {
        let parent_id = self.nodes[node_id].parent;
        let grandparent_id = self.nodes[parent_id].parent;
        self.nodes[parent_id].color = Color::Black;
        self.nodes[grandparent_id].color = Color::Red;
        if node_id == self.nodes[parent_id].left && parent_id == self.nodes[grandparent_id].left {
            self.rotate_right(parent_id);
        } else {
            self.rotate_left(parent_id);
        }
    }

    fn twin(&mut self, node_id: usize, parent_id: usize) -> usize {
        let parent_node = &self.nodes[parent_id];
        if node_id == parent_node.left {
            parent_node.right
        } else {
            parent_node.left
        }
    }

    fn rotate_left(&mut self, node_id: usize) {
        unsafe {
            let nodes_ptr = self.nodes.as_mut_ptr();
            let node = &mut *nodes_ptr.add(node_id);
            let parent_id = node.parent;
            let parent = &mut *nodes_ptr.add(parent_id);
            let grandparent_id = parent.parent;

            parent.right = node.left;
            if node.left != INVALID_IND {
                let left_child = &mut *nodes_ptr.add(node.left);
                left_child.parent = parent_id;
            }
            node.parent = grandparent_id;
            if grandparent_id != INVALID_IND {
                let grandparent = &mut *nodes_ptr.add(grandparent_id);
                if parent_id == grandparent.left {
                    grandparent.left = node_id;
                } else {
                    grandparent.right = node_id;
                }
            } else {
                self.root = node_id;
            }

            node.left = parent_id;
            parent.parent = node_id;
        }
    }

    fn rotate_right(&mut self, node_id: usize) {
        unsafe {
            let nodes_ptr = self.nodes.as_mut_ptr();
            let node = &mut *nodes_ptr.add(node_id);
            let parent_id = node.parent;
            let parent = &mut *nodes_ptr.add(parent_id);
            let grandparent_id = parent.parent;

            parent.left = node.right;
            if node.right != INVALID_IND {
                let right_child = &mut *nodes_ptr.add(node.right);
                right_child.parent = parent_id;
            }
            node.parent = grandparent_id;
            if grandparent_id != INVALID_IND {
                let grandparent = &mut *nodes_ptr.add(grandparent_id);
                if parent_id == grandparent.left {
                    grandparent.left = node_id;
                } else {
                    grandparent.right = node_id;
                }
            } else {
                self.root = node_id;
            }

            node.right = parent_id;
            parent.parent = node_id;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{alloc::Allocator, random};

    use crate::{INVALID_IND, utils::rbtree::Color};

    use super::RBTree;

    fn validate_tree<V: PartialOrd, A: Allocator + Copy>(tree: &RBTree<V, A>) -> bool {
        if tree.root == INVALID_IND {
            return true;
        }
        let mut stack = Vec::new();
        stack.push(tree.root);
        while let Some(node_id) = stack.pop() {
            let node = &tree.nodes[node_id];
            if node.color == super::Color::Red {
                if node.left != INVALID_IND && tree.nodes[node.left].color == super::Color::Red {
                    return false; // Red node with red child
                }
                if node.right != INVALID_IND && tree.nodes[node.right].color == super::Color::Red {
                    return false; // Red node with red child
                }
            }
            if node.left != INVALID_IND {
                stack.push(node.left);
            }
            if node.right != INVALID_IND {
                stack.push(node.right);
            }
        }
        black_height(tree, tree.root);
        true
    }

    fn black_height<V: PartialOrd, A: Allocator + Copy>(
        tree: &RBTree<V, A>,
        node_id: usize,
    ) -> usize {
        if node_id == INVALID_IND {
            return 1; // Null nodes contribute one black height
        }
        let left_height = black_height(tree, tree.nodes[node_id].left);
        let right_height = black_height(tree, tree.nodes[node_id].right);
        if tree.nodes[node_id].color == Color::Black && left_height != right_height {
            panic!("Black heights do not match at node {}", node_id);
        }
        if tree.nodes[node_id].color == super::Color::Black {
            left_height + 1
        } else {
            left_height
        }
    }

    #[test]
    fn test_insert() {
        let mut tree = RBTree::new(std::alloc::Global, false);
        // randomly generate 100 integers and insert them into the tree
        let values: Vec<i32> = (0..100).map(|_| random::random::<i32>() % 100).collect();
        for val in values {
            tree.insert(val);
            assert!(
                validate_tree(&tree),
                "Tree is not valid after inserting {}",
                val
            );
        }
    }

    #[test]
    fn test_range() {
        let mut tree = RBTree::new(std::alloc::Global, false);
        let values = [5, 3, 7, 2, 4, 6, 8, 1, 9];

        // Insert all values
        for &val in &values {
            tree.insert(val);
        }

        // Test full range iteration
        let items: Vec<i32> = tree.inorder_iter().copied().collect();
        assert_eq!(items, [1, 2, 3, 4, 5, 6, 7, 8, 9]);

        // Test range from 3..=7
        let range_items: Vec<i32> = tree.range(3..=7).copied().collect();
        assert_eq!(range_items, [3, 4, 5, 6, 7]);

        // Test range from 3..7 (exclusive upper bound)
        let range_items: Vec<i32> = tree.range(3..7).copied().collect();
        assert_eq!(range_items, [3, 4, 5, 6]);

        // Test unbounded ranges
        let range_items: Vec<i32> = tree.range(5..).copied().collect();
        assert_eq!(range_items, [5, 6, 7, 8, 9]);

        let range_items: Vec<i32> = tree.range(..5).copied().collect();
        assert_eq!(range_items, [1, 2, 3, 4]);
    }

    #[test]
    fn test_edge_cases() {
        // Test empty tree
        let empty_tree: RBTree<i32, _> = RBTree::new(std::alloc::Global, false);
        let items: Vec<i32> = empty_tree.inorder_iter().copied().collect();
        assert_eq!(items, []);

        let range_items: Vec<i32> = empty_tree.range(1..10).copied().collect();
        assert_eq!(range_items, []);

        // Test single node tree
        let mut single_tree = RBTree::new(std::alloc::Global, false);
        single_tree.insert(5);

        let items: Vec<i32> = single_tree.inorder_iter().copied().collect();
        assert_eq!(items, [5]);

        let range_items: Vec<i32> = single_tree.range(5..=5).copied().collect();
        assert_eq!(range_items, [5]);

        let range_items: Vec<i32> = single_tree.range(1..5).copied().collect();
        assert_eq!(range_items, []);

        let range_items: Vec<i32> = single_tree.range(6..10).copied().collect();
        assert_eq!(range_items, []);

        // Test larger tree with comprehensive ranges
        let mut tree = RBTree::new(std::alloc::Global, false);
        let values = [10, 5, 15, 3, 7, 12, 18, 1, 4, 6, 8, 11, 13, 16, 20];
        for &val in &values {
            tree.insert(val);
        }

        // Test ranges with no matches
        let range_items: Vec<i32> = tree.range(0..1).copied().collect();
        assert_eq!(range_items, []);

        let range_items: Vec<i32> = tree.range(21..30).copied().collect();
        assert_eq!(range_items, []);

        // Test ranges at exact boundaries
        let range_items: Vec<i32> = tree.range(1..=1).copied().collect();
        assert_eq!(range_items, [1]);

        let range_items: Vec<i32> = tree.range(20..=20).copied().collect();
        assert_eq!(range_items, [20]);

        // Test exclusive vs inclusive bounds
        let mut range_items: Vec<i32> = tree.range(5..10).copied().collect();
        range_items.sort();
        assert_eq!(range_items, [5, 6, 7, 8]);

        let mut range_items: Vec<i32> = tree.range(5..=10).copied().collect();
        range_items.sort();
        assert_eq!(range_items, [5, 6, 7, 8, 10]);

        // Test that iteration order is correct (in-order traversal)
        let items: Vec<i32> = tree.inorder_iter().copied().collect();
        let mut sorted_values = values.to_vec();
        sorted_values.sort();
        assert_eq!(items, sorted_values);

        // Test range with all elements
        let mut range_items: Vec<i32> = tree.range(..).copied().collect();
        range_items.sort();
        assert_eq!(range_items, sorted_values);
    }

    #[test]
    fn test_range_mut() {
        let mut tree = RBTree::new(std::alloc::Global, false);
        let values = [5, 3, 7, 2, 4, 6, 8, 1, 9];

        // Insert all values
        for &val in &values {
            tree.insert(val);
        }

        // Test basic range_mut functionality
        let mut count = 0;
        let mut found_values = Vec::new();
        for value in tree.range_mut(3..=7) {
            found_values.push(*value);
            *value += 100; // Add 100 to each value
            count += 1;
        }

        assert_eq!(count, 5); // Should find 5 values
        assert_eq!(found_values, [3, 4, 5, 6, 7]);

        // Test empty range
        let mut tree2 = RBTree::new(std::alloc::Global, false);
        tree2.insert(5);

        let mut empty_count = 0;
        for _value in tree2.range_mut(10..20) {
            empty_count += 1;
        }
        assert_eq!(empty_count, 0); // No values in range 10..20

        // Test single element modification
        let mut single_tree = RBTree::new(std::alloc::Global, false);
        single_tree.insert(42);

        for value in single_tree.range_mut(42..=42) {
            *value = 100;
        }

        let result: Vec<i32> = single_tree.inorder_iter().copied().collect();
        assert_eq!(result, [100]);
    }
}
