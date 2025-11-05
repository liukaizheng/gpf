#pragma once

#include <cmath>
#include <algorithm>
#include <numeric>

namespace gpf {
namespace math {

template<typename T>
inline void sub(const T* a, const T* b, T* c, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
}

inline double norm(const double* a, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum += a[i] * a[i];
    }
    return std::sqrt(sum);
}

template<typename T>
inline T dot(const T* a, const T* b, size_t n) {
    T ret = a[0] * b[0];
    for (size_t i = 1; i < n; ++i) {
        ret = ret + a[i] * b[i];
    }
    return ret;
}

template<typename T>
inline void cross(const T* a, const T* b, T* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

} // namespace math
} // namespace gpf
