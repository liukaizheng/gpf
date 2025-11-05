#pragma once

#include <vector>
#include <tuple>

namespace gpf {
namespace triangle {

struct DelaunayTriangulation {
    std::vector<double> points;
    std::vector<size_t> triangles;
    
    DelaunayTriangulation(const std::vector<double>& pts, const std::vector<std::pair<size_t, size_t>>& constraints = {});
    
    void add_constraint(size_t i, size_t j);
    void triangulate();
};

struct Tetrahedralization {
    std::vector<double> points;
    std::vector<size_t> tetrahedra;
    
    Tetrahedralization(const std::vector<double>& pts);
    
    void tetrahedralize();
};

std::tuple<std::vector<double>, std::vector<size_t>> 
triangulate_2d(const std::vector<double>& points, 
               const std::vector<std::pair<size_t, size_t>>& constraints = {});

std::tuple<std::vector<double>, std::vector<size_t>> 
tetrahedralize_3d(const std::vector<double>& points);

} // namespace triangle
} // namespace gpf
