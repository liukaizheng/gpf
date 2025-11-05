#pragma once

#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>

namespace gpf {
namespace predicates {

enum class Orientation {
    Positive,
    Negative,
    Zero,
    Undefined
};

inline Orientation double_to_sign(double x) {
    if (x > 0.0) return Orientation::Positive;
    if (x < 0.0) return Orientation::Negative;
    return Orientation::Zero;
}

inline Orientation sign_reverse(Orientation ori) {
    if (ori == Orientation::Positive) return Orientation::Negative;
    if (ori == Orientation::Negative) return Orientation::Positive;
    return ori;
}

inline bool sign_reversed(Orientation ori1, Orientation ori2) {
    return (ori1 == Orientation::Positive && ori2 == Orientation::Negative) ||
           (ori2 == Orientation::Positive && ori1 == Orientation::Negative);
}

inline int get_exponent(double x) {
    if (x == 0.0) return 0;
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof(double));
    return static_cast<int>((bits >> 52) & 0x7FF) - 1023;
}

double orient2d(const double* pa, const double* pb, const double* pc);

double orient3d(const double* pa, const double* pb, const double* pc, const double* pd);

double incircle(const double* pa, const double* pb, const double* pc, const double* pd);

double insphere(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe);

inline bool same_point(const double* p, const double* q) {
    return p[0] == q[0] && p[1] == q[1] && p[2] == q[2];
}

bool point_in_inner_triangle(const double* p, const double* v1, const double* v2, const double* v3);

bool inner_segment_cross_inner_triangle(const double* u1, const double* u2,
                                        const double* v1, const double* v2, const double* v3);

bool inner_segments_cross(const double* u1, const double* u2,
                          const double* v1, const double* v2);

bool point_in_inner_segment(const double* p, const double* v1, const double* v2);

size_t max_comp_in_tri_normal(const double* ov1, const double* ov2, const double* ov3);

bool same_half_plane(const double* p, const double* q, const double* v1, const double* v2);

bool mis_alignment(const double* p, const double* q, const double* r);

} // namespace predicates
} // namespace gpf
