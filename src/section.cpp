#include "section.h"

#include "force.h"
#include "polygon.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <random>

#include <png.h>

#ifdef FRAPS_USE_EIGEN
template<typename VectorA, typename VectorB>
Scalar dot(VectorA a, VectorB b) {
    return a.dot(b);
}

template<typename VectorA, typename VectorB>
Scalar dist(VectorA a, VectorB b) {
    return (a - b).norm();
}

template<typename Vector>
void normalize(Vector a) {
    a.normalize();
}
#endif

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

#ifndef NDEBUG
    std::cout << "dpde " << dpde[0] << " " << dpde[1] << std::endl;
    std::cout << "p " << p[0] << " " << p[1] << std::endl;
#endif

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

#ifndef NDEBUG
    std::cout << "section found" << std::endl;
    std::cout << "v0" << ret.j0.p[0] << " " << ret.j0.p[1] << std::endl;
    std::cout << "v1" << ret.j1.p[0] << " " << ret.j1.p[1] << std::endl;
#endif

    normalize(ret.normal);
    // let normal points to s0
    if (ret.j0.idx_v0 < ret.j1.idx_v0) {
        if (dot(shape[i0.idx_v1] - shape[i0.idx_v0], ret.normal) > 0) {
            ret.normal = -ret.normal;
        }
    } else {
        if (dot(shape[i1.idx_v1] - shape[i1.idx_v0], ret.normal) > 0) {
            ret.normal = -ret.normal;
        }
    }

    return true;
}

std::vector<Section> find_cross_sections(const std::vector<Vector>& shape) {
    std::vector<Section> ret;
    const int num_directions = 12;
    const int num_offsets = 15;
    const double left = -1.0;
    const double right = 1.0;

    const double length_threshold = 0.01;

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
#ifndef NDEBUG
            std::cout << std::endl << normal[0] << " " << normal[1] << " " << offset << std::endl;
            for (auto const & itsc : intersections) {
                std::cout << itsc.p[0] << " " << itsc.p[1] << std::endl;
            }
#endif

            assert(intersections.size() % 2 == 0);

            std::sort(intersections.begin(), intersections.end(), [](const Intersection& a, const Intersection& b) {
                if (std::abs(a.p[0] - b.p[0]) < cross_section_align_tol) { return a.p[1] > b.p[1]; }
                else { return a.p[0] > b.p[0]; }
            });

            for (unsigned int k = 0; k < intersections.size(); k += 2) {
                Section s;
                if (construct_cross_section(intersections[k], intersections[k + 1], shape, s)) {
                    if (s.length > length_threshold)
                        ret.push_back(s);
                }
            }
        }
    }

    return ret;
}
