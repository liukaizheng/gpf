#include "gpf/math.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {4.0, 5.0, 6.0};
    double c[3];
    
    gpf::math::sub(a, b, c, 3);
    assert(c[0] == -3.0 && c[1] == -3.0 && c[2] == -3.0);
    
    double norm_val = gpf::math::norm(a, 3);
    assert(std::abs(norm_val - std::sqrt(14.0)) < 1e-10);
    std::cout << "Norm: " << norm_val << "\n";
    
    double dot_val = gpf::math::dot(a, b, 3);
    assert(dot_val == 32.0);
    std::cout << "Dot product: " << dot_val << "\n";
    
    double cross[3];
    gpf::math::cross(a, b, cross);
    assert(cross[0] == -3.0 && cross[1] == 6.0 && cross[2] == -3.0);
    
    std::cout << "All math tests passed!\n";
    return 0;
}
