#include <gpf/gpf.hpp>
#include <iostream>
#include <vector>

int main() {
    std::cout << "=== GPF C++ Library Demo ===\n\n";
    
    // 1. Math operations
    std::cout << "1. Math Operations\n";
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {4.0, 5.0, 6.0};
    double c[3];
    
    gpf::math::sub(a, b, c, 3);
    std::cout << "   Subtraction: [" << c[0] << ", " << c[1] << ", " << c[2] << "]\n";
    
    double dot = gpf::math::dot(a, b, 3);
    std::cout << "   Dot product: " << dot << "\n";
    
    double cross[3];
    gpf::math::cross(a, b, cross);
    std::cout << "   Cross product: [" << cross[0] << ", " << cross[1] << ", " << cross[2] << "]\n";
    std::cout << std::endl;
    
    // 2. Disjoint Set
    std::cout << "2. Disjoint Set (Union-Find)\n";
    gpf::DisjointSet ds(6);
    std::cout << "   Initial groups: " << ds.n_groups << "\n";
    
    ds.merge(0, 1);
    ds.merge(2, 3);
    ds.merge(4, 5);
    std::cout << "   After merging pairs: " << ds.n_groups << " groups\n";
    
    ds.merge(0, 4);
    std::cout << "   After connecting components: " << ds.n_groups << " groups\n";
    std::cout << std::endl;
    
    // 3. Graph Cut
    std::cout << "3. Graph Cut (Max Flow)\n";
    std::vector<double> source_cap = {10.0, 8.0, 12.0, 15.0};
    std::vector<double> sink_cap = {7.0, 9.0, 10.0, 11.0};
    
    gpf::GraphCut gc(source_cap, sink_cap);
    gc.add_edge(0, 1, 5.0, 5.0);
    gc.add_edge(1, 2, 8.0, 8.0);
    gc.add_edge(2, 3, 6.0, 6.0);
    gc.add_edge(0, 2, 10.0, 10.0);
    
    double flow = gc.max_flow();
    std::cout << "   Maximum flow: " << flow << "\n";
    std::cout << "   Nodes in source set: ";
    for (size_t i = 0; i < gc.is_sink.size(); ++i) {
        if (!gc.is_sink[i]) std::cout << i << " ";
    }
    std::cout << "\n";
    std::cout << "   Nodes in sink set: ";
    for (size_t i = 0; i < gc.is_sink.size(); ++i) {
        if (gc.is_sink[i]) std::cout << i << " ";
    }
    std::cout << "\n";
    std::cout << std::endl;
    
    // 4. Geometric Predicates
    std::cout << "4. Geometric Predicates\n";
    
    // 2D orientation
    double p1[] = {0.0, 0.0};
    double p2[] = {1.0, 0.0};
    double p3[] = {0.5, 0.5};
    
    double orient = gpf::predicates::orient2d(p1, p2, p3);
    std::cout << "   Orient2D(p1, p2, p3): " << orient << " (";
    if (orient > 0) std::cout << "counter-clockwise)\n";
    else if (orient < 0) std::cout << "clockwise)\n";
    else std::cout << "collinear)\n";
    
    // 3D orientation
    double q1[] = {0.0, 0.0, 0.0};
    double q2[] = {1.0, 0.0, 0.0};
    double q3[] = {0.0, 1.0, 0.0};
    double q4[] = {0.0, 0.0, 1.0};
    
    double orient3 = gpf::predicates::orient3d(q1, q2, q3, q4);
    std::cout << "   Orient3D(q1, q2, q3, q4): " << orient3 << "\n";
    
    // Incircle test
    double r1[] = {0.0, 0.0};
    double r2[] = {1.0, 0.0};
    double r3[] = {0.0, 1.0};
    double r4[] = {0.5, 0.5};
    
    double incircle = gpf::predicates::incircle(r1, r2, r3, r4);
    std::cout << "   InCircle test: " << incircle << " (";
    if (incircle > 0) std::cout << "inside)\n";
    else if (incircle < 0) std::cout << "outside)\n";
    else std::cout << "on circle)\n";
    
    std::cout << "\n=== Demo Complete ===\n";
    
    return 0;
}
