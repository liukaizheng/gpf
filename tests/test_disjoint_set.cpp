#include "gpf/disjoint_set.hpp"
#include <iostream>
#include <cassert>

int main() {
    gpf::DisjointSet ds(5);
    assert(ds.n_groups == 5);
    
    ds.merge(0, 1);
    assert(ds.n_groups == 4);
    
    ds.merge(2, 3);
    assert(ds.n_groups == 3);
    
    ds.merge(0, 2);
    assert(ds.n_groups == 2);
    
    auto groups = ds.output();
    assert(groups.size() == 2);
    
    std::cout << "All disjoint set tests passed!\n";
    return 0;
}
