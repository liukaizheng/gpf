#include "gpf/predicates.hpp"
#include <cmath>
#include <algorithm>

namespace gpf {
namespace predicates {

double orient2d(const double* pa, const double* pb, const double* pc) {
    double acx = pa[0] - pc[0];
    double bcx = pb[0] - pc[0];
    double acy = pa[1] - pc[1];
    double bcy = pb[1] - pc[1];
    return acx * bcy - acy * bcx;
}

double orient3d(const double* pa, const double* pb, const double* pc, const double* pd) {
    double adx = pa[0] - pd[0];
    double bdx = pb[0] - pd[0];
    double cdx = pc[0] - pd[0];
    double ady = pa[1] - pd[1];
    double bdy = pb[1] - pd[1];
    double cdy = pc[1] - pd[1];
    double adz = pa[2] - pd[2];
    double bdz = pb[2] - pd[2];
    double cdz = pc[2] - pd[2];
    
    return adx * (bdy * cdz - bdz * cdy)
         + bdx * (cdy * adz - cdz * ady)
         + cdx * (ady * bdz - adz * bdy);
}

double incircle(const double* pa, const double* pb, const double* pc, const double* pd) {
    double adx = pa[0] - pd[0];
    double ady = pa[1] - pd[1];
    double bdx = pb[0] - pd[0];
    double bdy = pb[1] - pd[1];
    double cdx = pc[0] - pd[0];
    double cdy = pc[1] - pd[1];
    
    double abdet = adx * bdy - bdx * ady;
    double bcdet = bdx * cdy - cdx * bdy;
    double cadet = cdx * ady - adx * cdy;
    double alift = adx * adx + ady * ady;
    double blift = bdx * bdx + bdy * bdy;
    double clift = cdx * cdx + cdy * cdy;
    
    return alift * bcdet + blift * cadet + clift * abdet;
}

double insphere(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe) {
    double aex = pa[0] - pe[0];
    double bex = pb[0] - pe[0];
    double cex = pc[0] - pe[0];
    double dex = pd[0] - pe[0];
    double aey = pa[1] - pe[1];
    double bey = pb[1] - pe[1];
    double cey = pc[1] - pe[1];
    double dey = pd[1] - pe[1];
    double aez = pa[2] - pe[2];
    double bez = pb[2] - pe[2];
    double cez = pc[2] - pe[2];
    double dez = pd[2] - pe[2];
    
    double ab = aex * bey - bex * aey;
    double bc = bex * cey - cex * bey;
    double cd = cex * dey - dex * cey;
    double da = dex * aey - aex * dey;
    
    double ac = aex * cey - cex * aey;
    double bd = bex * dey - dex * bey;
    
    double abc = aez * bc - bez * ac + cez * ab;
    double bcd = bez * cd - cez * bd + dez * bc;
    double cda = cez * da + dez * ac + aez * cd;
    double dab = dez * ab + aez * bd + bez * da;
    
    double alift = aex * aex + aey * aey + aez * aez;
    double blift = bex * bex + bey * bey + bez * bez;
    double clift = cex * cex + cey * cey + cez * cez;
    double dlift = dex * dex + dey * dey + dez * dez;
    
    return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
}

bool point_in_inner_triangle(const double* p, const double* v1, const double* v2, const double* v3) {
    auto o1 = double_to_sign(orient2d(p, v2, v3));
    auto o2 = double_to_sign(orient2d(v1, v2, v3));
    if (o1 != o2) return false;
    
    double p_yz[2] = {p[1], p[2]};
    double v1_yz[2] = {v1[1], v1[2]};
    double v2_yz[2] = {v2[1], v2[2]};
    double v3_yz[2] = {v3[1], v3[2]};
    
    o1 = double_to_sign(orient2d(p_yz, v2_yz, v3_yz));
    o2 = double_to_sign(orient2d(v1_yz, v2_yz, v3_yz));
    if (o1 != o2) return false;
    
    double p_xz[2] = {p[0], p[2]};
    double v1_xz[2] = {v1[0], v1[2]};
    double v2_xz[2] = {v2[0], v2[2]};
    double v3_xz[2] = {v3[0], v3[2]};
    
    o1 = double_to_sign(orient2d(p_xz, v2_xz, v3_xz));
    o2 = double_to_sign(orient2d(v1_xz, v2_xz, v3_xz));
    if (o1 != o2) return false;
    
    return true;
}

bool inner_segment_cross_inner_triangle(const double* u1, const double* u2,
                                        const double* v1, const double* v2, const double* v3) {
    double bound = std::min(u1[0], u2[0]);
    if (v1[0] <= bound && v2[0] <= bound && v3[0] <= bound) return false;
    
    bound = std::max(u1[0], u2[0]);
    if (v1[0] >= bound && v2[0] >= bound && v3[0] >= bound) return false;
    
    bound = std::min(u1[1], u2[1]);
    if (v1[1] <= bound && v2[1] <= bound && v3[1] <= bound) return false;
    
    bound = std::max(u1[1], u2[1]);
    if (v1[1] >= bound && v2[1] >= bound && v3[1] >= bound) return false;
    
    bound = std::min(u1[2], u2[2]);
    if (v1[2] <= bound && v2[2] <= bound && v3[2] <= bound) return false;
    
    bound = std::max(u1[2], u2[2]);
    if (v1[2] >= bound && v2[2] >= bound && v3[2] >= bound) return false;
    
    auto orient_u1_tri = double_to_sign(orient3d(u1, v1, v2, v3));
    auto orient_u2_tri = double_to_sign(orient3d(u2, v1, v2, v3));
    
    if (orient_u1_tri == Orientation::Zero || orient_u2_tri == Orientation::Zero) return false;
    if (orient_u1_tri == orient_u2_tri) return false;
    
    auto orient_u_v1v2 = double_to_sign(orient3d(u1, u2, v1, v2));
    auto orient_u_v2v3 = double_to_sign(orient3d(u1, u2, v2, v3));
    
    if (orient_u_v1v2 == Orientation::Zero || orient_u_v2v3 == Orientation::Zero) return false;
    if (orient_u_v1v2 != orient_u_v2v3) return false;
    
    auto orient_u_v3v1 = double_to_sign(orient3d(u1, u2, v3, v1));
    if (orient_u_v3v1 == Orientation::Zero) return false;
    if (orient_u_v3v1 != orient_u_v2v3) return false;
    
    return true;
}

bool same_half_plane(const double* p, const double* q, const double* v1, const double* v2) {
    if (double_to_sign(orient2d(p, v1, v2)) != double_to_sign(orient2d(q, v1, v2))) return false;
    
    double p_yz[2] = {p[1], p[2]};
    double q_yz[2] = {q[1], q[2]};
    double v1_yz[2] = {v1[1], v1[2]};
    double v2_yz[2] = {v2[1], v2[2]};
    
    if (double_to_sign(orient2d(p_yz, v1_yz, v2_yz)) != double_to_sign(orient2d(q_yz, v1_yz, v2_yz))) return false;
    
    double p_xz[2] = {p[0], p[2]};
    double q_xz[2] = {q[0], q[2]};
    double v1_xz[2] = {v1[0], v1[2]};
    double v2_xz[2] = {v2[0], v2[2]};
    
    return double_to_sign(orient2d(p_xz, v1_xz, v2_xz)) == double_to_sign(orient2d(q_xz, v1_xz, v2_xz));
}

bool mis_alignment(const double* p, const double* q, const double* r) {
    if (orient2d(p, q, r) != 0.0) return true;
    
    double p_yz[2] = {p[1], p[2]};
    double q_yz[2] = {q[1], q[2]};
    double r_yz[2] = {r[1], r[2]};
    if (orient2d(p_yz, q_yz, r_yz) != 0.0) return true;
    
    double p_xz[2] = {p[0], p[2]};
    double q_xz[2] = {q[0], q[2]};
    double r_xz[2] = {r[0], r[2]};
    return orient2d(p_xz, q_xz, r_xz) != 0.0;
}

bool inner_segments_cross(const double* u1, const double* u2, const double* v1, const double* v2) {
    if (orient3d(u1, u2, v1, v2) != 0.0) return false;
    
    if (same_half_plane(u1, u2, v1, v2) || same_half_plane(v1, v2, u1, u2)) return false;
    
    if (!mis_alignment(u1, v1, v2)) return false;
    if (!mis_alignment(u2, v1, v2)) return false;
    if (!mis_alignment(v1, u1, u2)) return false;
    if (!mis_alignment(v2, u1, u2)) return false;
    
    if (orient2d(u1, u2, v1) != 0.0) return true;
    if (orient2d(v1, v2, u2) != 0.0) return true;
    
    double u1_yz[2] = {u1[1], u1[2]};
    double u2_yz[2] = {u2[1], u2[2]};
    double v1_yz[2] = {v1[1], v1[2]};
    double v2_yz[2] = {v2[1], v2[2]};
    
    if (orient2d(u1_yz, u2_yz, v1_yz) != 0.0) return true;
    if (orient2d(v1_yz, v2_yz, u2_yz) != 0.0) return true;
    
    double u1_xz[2] = {u1[0], u1[2]};
    double u2_xz[2] = {u2[0], u2[2]};
    double v1_xz[2] = {v1[0], v1[2]};
    double v2_xz[2] = {v2[0], v2[2]};
    
    if (orient2d(u1_xz, u2_xz, v1_xz) != 0.0) return true;
    if (orient2d(v1_xz, v2_xz, u2_xz) != 0.0) return true;
    
    return false;
}

bool point_in_inner_segment(const double* p, const double* v1, const double* v2) {
    if (mis_alignment(p, v1, v2)) return false;
    
    return ((v1[0] < v2[0] && v1[0] < p[0] && p[0] < v2[0]) ||
            (v1[0] > v2[0] && v1[0] > p[0] && p[0] > v2[0]) ||
            (v1[1] < v2[1] && v1[1] < p[1] && p[1] < v2[1]) ||
            (v1[1] > v2[1] && v1[1] > p[1] && p[1] > v2[1]) ||
            (v1[2] < v2[2] && v1[2] < p[2] && p[2] < v2[2]) ||
            (v1[2] > v2[2] && v1[2] > p[2] && p[2] > v2[2]));
}

size_t max_comp_in_tri_normal(const double* ov1, const double* ov2, const double* ov3) {
    double v3x = ov3[0] - ov2[0];
    double v3y = ov3[1] - ov2[1];
    double v3z = ov3[2] - ov2[2];
    double v2x = ov2[0] - ov1[0];
    double v2y = ov2[1] - ov1[1];
    double v2z = ov2[2] - ov1[2];
    
    double nvx = v2y * v3z - v2z * v3y;
    double nvy = v3x * v2z - v3z * v2x;
    double nvz = v2x * v3y - v2y * v3x;
    
    double abs_nvx = std::abs(nvx);
    double abs_nvy = std::abs(nvy);
    double abs_nvz = std::abs(nvz);
    
    if (abs_nvx >= abs_nvy && abs_nvx >= abs_nvz) return 0;
    if (abs_nvy >= abs_nvz) return 1;
    return 2;
}

} // namespace predicates
} // namespace gpf
