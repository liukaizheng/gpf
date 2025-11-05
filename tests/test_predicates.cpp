#include "gpf/predicates.hpp"
#include <iostream>
#include <cassert>

int main() {
    double pa[] = {0.0, 0.0};
    double pb[] = {1.0, 0.0};
    double pc[] = {0.0, 1.0};
    
    double o2d = gpf::predicates::orient2d(pa, pb, pc);
    assert(o2d > 0.0);
    std::cout << "Orient2D result: " << o2d << "\n";
    
    double p3a[] = {0.0, 0.0, 0.0};
    double p3b[] = {1.0, 0.0, 0.0};
    double p3c[] = {0.0, 1.0, 0.0};
    double p3d[] = {0.0, 0.0, 1.0};
    
    double o3d = gpf::predicates::orient3d(p3a, p3b, p3c, p3d);
    assert(o3d != 0.0);
    std::cout << "Orient3D result: " << o3d << "\n";
    
    double ic = gpf::predicates::incircle(pa, pb, pc, pb);
    std::cout << "Incircle result: " << ic << "\n";
    
    std::cout << "All predicate tests passed!\n";
    return 0;
}
