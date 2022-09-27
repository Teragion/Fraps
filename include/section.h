#ifndef SECTION_H
#define SECTION_H

#include "force.h"
#include "polygon.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <random>

#include <png.h>

constexpr double cross_section_align_tol = 0.000001;

using Scalar = double;
using Vector = Vec2d;

struct Intersection {
    Vector p;
    int idx_v0 = -1;
    int idx_v1 = -1;
};

struct Section {
    Intersection j0;
    Intersection j1;
    Vector p;
    Vector normal;
    Scalar length = 0.0;

    // always true that dot(c1, normal) > 0
    Vector c0;
    Vector c1;

    Scalar I0;
    Scalar I0_inv{};
    Scalar I1;
    Scalar I1_inv{};

    Scalar w0;
    Scalar w0_inv{};
    Scalar w1;
    Scalar w1_inv{};
};

/**
 * @brief finds intersection between a line (described by normal and offset) and
 * a line segment
 *
 * @param normal
 * @param offset
 * @param v1
 * @param v2
 * @return Vector
 */
bool intersection_normal_segment(const Vector& normal, const double& offset,
                                 const Vector& v0, const Vector& v1, Vector& ret);

bool construct_cross_section(const Intersection& i0, const Intersection& i1, const Polygon<double, Vector >& shape, Section& ret);

/**
 * @brief constructs a vector of sections from a rigidpoly
 * @param shape
 * @return
 */
std::vector<Section> find_cross_sections(const std::vector<Vector>& shape);


inline double moment_of_inertia(double l) {
    double r = l / 2;
    return 2 * (r * r * r / 3);
}

bool intersection_normal_segment(const Vector& normal, const double& offset,
                                 const Vector& v0, const Vector& v1, Vector& ret) {
    bool p0 = dot(v0, normal) < offset;
    bool p1 = dot(v1, normal) < offset;

    if (p0 == p1) { return false; }

    Vector dir = v0 - v1;
    Vector normal2 = {dir[1], -dir[0]};
    double offset2 = dot(normal2, v0);

    ret = {(offset * normal2[1] - normal[1] * offset2) / (normal[0] * normal2[1] - normal[1] * normal2[0]),
           (offset * normal2[0] - normal[0] * offset2) / (normal[1] * normal2[0] - normal[0] * normal2[1])};

    return true;
}

inline Vector find_edge_normal(unsigned int idx_v0, unsigned int idx_v1, const std::vector<Vector>& shape) {
    Vector ret{-(shape[idx_v0][1] - shape[idx_v1][1]), shape[idx_v0][0] - shape[idx_v1][0]};
    normalize(ret);
    return ret;
}

Vector find_advance_pd(const Vector& n, const Vector& n0) {
    double det_inv = 1.0 / (n0[0] * n[1] - n0[1] * n[0]);
    return {det_inv * (-n0[1]), det_inv * n0[0]};
}

bool construct_cross_section(const Intersection& i0, const Intersection& i1, const std::vector<Vector>& shape, Section& ret) {
    // find dpde
    Vector n{-(i0.p[1] - i1.p[1]), i0.p[0] - i1.p[0]}; // normal of line of raw intersection
    Vector n0 = find_edge_normal(i0.idx_v0, i0.idx_v1, shape);
    Vector n1 = find_edge_normal(i1.idx_v0, i1.idx_v1, shape);

    if (dot(n, n0) < 0) { n0 = -n0; }
    if (dot(n, n1) < 0) { n1 = -n1; }

    if (std::abs(n0[0] * n[1] - n0[1] * n[0]) < cross_section_align_tol) { return false; }
    Vector dqde = find_advance_pd(n, n0);

    if (std::abs(n1[0] * n[1] - n1[1] * n[0]) < cross_section_align_tol) { return false; }
    Vector drde = find_advance_pd(n, n1);

    Vector dpde = (dqde + drde) * 0.5;

    Vector p = (i0.p + i1.p) * 0.5;

    std::cout << "dpde " << dpde[0] << " " << dpde[1] << std::endl;
    std::cout << "p " << p[0] << " " << p[1] << std::endl;

    ret.normal = dpde;

    // reconstruct intersections
    ret.j0.idx_v0 = i0.idx_v0;
    ret.j0.idx_v1 = i0.idx_v1;
    if (!intersection_normal_segment(dpde, dot(dpde, p), shape[i0.idx_v0], shape[i0.idx_v1], ret.j0.p)) { return false; }

    ret.j1.idx_v0 = i1.idx_v0;
    ret.j1.idx_v1 = i1.idx_v1;
    if (!intersection_normal_segment(dpde, dot(dpde, p), shape[i1.idx_v0], shape[i1.idx_v1], ret.j1.p)) { return false; }

    ret.length = dist(ret.j0.p, ret.j1.p);
    ret.p = p;

    std::cout << "section found" << std::endl;
    std::cout << "v0" << ret.j0.p[0] << " " << ret.j0.p[1] << std::endl;
    std::cout << "v1" << ret.j1.p[0] << " " << ret.j1.p[1] << std::endl;

    return true;
}

std::vector<Section> find_cross_sections(const std::vector<Vector>& shape) {
    std::vector<Section> ret;
    const int num_directions = 12;
    const int num_offsets = 1;
    const double left = -0.0;
    const double right = 0.0;

    for (int i = 0; i < num_directions; i++) {
        const double theta = i * (M_PI / static_cast<double>(num_directions));
        for (int j = 0; j < num_offsets; j++) {
            const double offset = lerp(left, right, j / static_cast<double>(num_offsets));

            const Vector normal{std::cos(theta), std::sin(theta)};

            std::vector<Intersection> intersections;

            for (int k = 0; k < shape.size(); k++) {
                Vector p;
                if (!intersection_normal_segment(normal, offset, shape[k],
                                                 shape[(k + 1) % shape.size()], p)) {
                    continue;
                }

                intersections.push_back({p, k, static_cast<int>((k + 1) % shape.size())});
            }
            std::cout << std::endl << normal[0] << " " << normal[1] << " " << offset << std::endl;
            for (auto const & itsc : intersections) {
                std::cout << itsc.p[0] << " " << itsc.p[1] << std::endl;
            }

            assert(intersections.size() % 2 == 0);

            std::sort(intersections.begin(), intersections.end(), [](const Intersection& a, const Intersection& b) {
                if (std::abs(a.p[0] - b.p[0]) < cross_section_align_tol) { return a.p[1] > b.p[1]; }
                else { return a.p[0] > b.p[0]; }
            });

            for (unsigned int k = 0; k < intersections.size(); k += 2) {
                Section s;
                if (construct_cross_section(intersections[k], intersections[k + 1], shape, s)) {
                    ret.push_back(s);
                }
            }
        }
    }

    return ret;
}

#endif  // SECTION_H
