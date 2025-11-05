#include "gpf/triangle.hpp"
#include "gpf/predicates.hpp"
#include <algorithm>
#include <stdexcept>

namespace gpf {
namespace triangle {

DelaunayTriangulation::DelaunayTriangulation(const std::vector<double>& pts, const std::vector<std::pair<size_t, size_t>>& constraints) 
    : points(pts) {
}

void DelaunayTriangulation::add_constraint(size_t i, size_t j) {
}

void DelaunayTriangulation::triangulate() {
}

Tetrahedralization::Tetrahedralization(const std::vector<double>& pts)
    : points(pts) {
}

void Tetrahedralization::tetrahedralize() {
}

std::tuple<std::vector<double>, std::vector<size_t>> 
triangulate_2d(const std::vector<double>& points, const std::vector<std::pair<size_t, size_t>>& constraints) {
    return {points, {}};
}

std::tuple<std::vector<double>, std::vector<size_t>> 
tetrahedralize_3d(const std::vector<double>& points) {
    return {points, {}};
}

} // namespace triangle
} // namespace gpf
