#pragma once

#include "gpf/math.hpp"
#include "gpf/disjoint_set.hpp"
#include "gpf/graphcut.hpp"
#include "gpf/predicates.hpp"
#include "gpf/triangle.hpp"
#include "gpf/mesh.hpp"
#include "gpf/polygonlization.hpp"

namespace gpf {

constexpr size_t INVALID_IND = std::numeric_limits<size_t>::max();

inline const double* point(const double* points, size_t tid) {
    return &points[tid * 3];
}

inline double face_area_2d(const std::vector<double>& points) {
    double area = 0.0;
    size_t n = points.size() / 2;
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        area += points[i * 2] * points[j * 2 + 1] - points[i * 2 + 1] * points[j * 2];
    }
    return area;
}

} // namespace gpf
