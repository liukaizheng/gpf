#include "gpf/polygonlization.hpp"
#include "gpf/predicates.hpp"
#include "gpf/disjoint_set.hpp"
#include <unordered_map>
#include <algorithm>

namespace gpf {
namespace polygonlization {

std::tuple<std::vector<double>, std::vector<size_t>> 
make_mesh_for_triangles(const std::vector<double>& points,
                        const std::vector<size_t>& triangles,
                        const std::vector<size_t>& tri_in_shells) {
    (void)tri_in_shells;
    return {points, triangles};
}

BSPComplex::BSPComplex(const std::vector<double>& points, const std::vector<size_t>& triangles)
    : points_(points), triangles_(triangles) {}

void BSPComplex::split_intersecting_faces() {
}

void BSPComplex::remove_duplicates() {
}

std::tuple<std::vector<double>, std::vector<size_t>> BSPComplex::extract_mesh() {
    return {points_, triangles_};
}

ConformingMesh::ConformingMesh(const std::vector<double>& points, const std::vector<size_t>& triangles)
    : points_(points), triangles_(triangles) {}

void ConformingMesh::make_conforming() {
}

std::tuple<std::vector<double>, std::vector<size_t>> ConformingMesh::get_mesh() {
    return {points_, triangles_};
}

} // namespace polygonlization
} // namespace gpf
