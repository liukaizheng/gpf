#pragma once

#include <vector>
#include <unordered_map>
#include <numeric>

namespace gpf {

class DisjointSet {
public:
    explicit DisjointSet(size_t n) 
        : n_groups(n), parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    void merge(size_t x, size_t y) {
        size_t root_x = find_set(x);
        size_t root_y = find_set(y);
        link(root_x, root_y);
    }

    std::unordered_map<size_t, std::vector<size_t>> output() {
        std::unordered_map<size_t, std::vector<size_t>> result;
        for (size_t i = 0; i < parent_.size(); ++i) {
            size_t root = find_set(i);
            result[root].push_back(i);
        }
        return result;
    }

    size_t n_groups;

private:
    std::vector<size_t> parent_;
    std::vector<size_t> rank_;

    void link(size_t x, size_t y) {
        if (x == y) {
            return;
        }

        if (rank_[x] > rank_[y]) {
            parent_[y] = x;
        } else {
            parent_[x] = y;
            if (rank_[x] == rank_[y]) {
                rank_[y]++;
            }
        }
        n_groups--;
    }

    size_t find_set(size_t x) {
        if (x != parent_[x]) {
            parent_[x] = find_set(parent_[x]);
        }
        return parent_[x];
    }
};

} // namespace gpf
