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

#include <Eigen/Dense>

constexpr double cross_section_align_tol = 0.000001;

using Scalar = double;
#ifndef FRAPS_USE_EIGEN
using Vector = Vec2d;
#else
using Vector = Eigen::Matrix<Scalar, 3, 1>;
#endif

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

    Scalar ratio = 0.0;

    int fixed_side = -1;
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

inline Vector find_edge_normal(unsigned int idx_v0, unsigned int idx_v1, const std::vector<Vector>& shape) {
    Vector ret{-(shape[idx_v0][1] - shape[idx_v1][1]), shape[idx_v0][0] - shape[idx_v1][0]};
#ifndef FRAPS_USE_EIGEN
    normalize(ret);
#else
    ret.normalize();
#endif
    return ret;
}

Vector find_advance_pd(const Vector& n, const Vector& n0);

bool construct_cross_section(const Intersection& i0, const Intersection& i1, const std::vector<Vector>& shape, Section& ret);

std::vector<Section> find_cross_sections(const std::vector<Vector>& shape);

#endif // SECTION_H
