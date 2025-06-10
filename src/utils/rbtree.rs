use std::{alloc::Allocator, collections::VecDeque};

use crate::INVALID_IND;

#[derive(PartialEq, Eq)]
enum Color {
    Red,
    Black,
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
        let uncle_id = self.uncle(node_id);
        if uncle_id != INVALID_IND && self.nodes[uncle_id].color == Color::Red {
            self.nodes[node_id].color = Color::Black;
            self.nodes[uncle_id].color = Color::Black;
            let grandparent_id = self.grandparent(node_id);
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

    fn grandparent(&mut self, node_id: usize) -> usize {
        let parent_id = self.nodes[node_id].parent;
        if parent_id == INVALID_IND {
            INVALID_IND
        } else {
            self.nodes[parent_id].parent
        }
    }

    fn uncle(&mut self, node_id: usize) -> usize {
        let grandparent_id = self.grandparent(node_id);
        if grandparent_id == INVALID_IND {
            INVALID_IND
        } else {
            let grandparent_node = &self.nodes[grandparent_id];
            if self.nodes[node_id].parent == grandparent_node.left {
                grandparent_node.right
            } else {
                grandparent_node.left
            }
        }
    }

    fn rotate_left(&mut self, node_id: usize) {
        unsafe {
            let nodes_ptr = self.nodes.as_mut_ptr();
            let node = &mut *nodes_ptr.add(node_id);
            let parent_id = node.parent;
            let parent = &mut *nodes_ptr.add(parent_id);
            let grandparent_id = parent.parent;
            let grandparent = &mut *nodes_ptr.add(grandparent_id);

            parent.right = node.left;
            if node.left != INVALID_IND {
                let left_child = &mut *nodes_ptr.add(node.left);
                left_child.parent = parent_id;
            }
            node.parent = grandparent_id;
            if parent_id == grandparent.left {
                grandparent.left = node_id;
            } else {
                grandparent.right = node_id;
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
            let grandparent = &mut *nodes_ptr.add(grandparent_id);

            parent.left = node.right;
            if node.right != INVALID_IND {
                let right_child = &mut *nodes_ptr.add(node.right);
                right_child.parent = parent_id;
            }
            node.parent = grandparent_id;
            if parent_id == grandparent.left {
                grandparent.left = node_id;
            } else {
                grandparent.right = node_id;
            }

            node.right = parent_id;
            parent.parent = node_id;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::alloc::Allocator;

    use crate::INVALID_IND;

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
        true
    }

    #[test]
    fn test_insert() {
        let mut tree = RBTree::new(std::alloc::Global, false);
        for v in [2, 3, 1, 0, 7, 9, 4, 8, 5, 6] {
            tree.insert(v);
            assert!(validate_tree(&tree), "Tree is not valid after inserting {}", v);
        }
        let a = 2;
    }
}
