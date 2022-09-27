#ifndef SAT_H
#define SAT_H

#include "collision_parameters.h"
#include "vec.h"

#include <complex>
#include <limits>
#include <vector>
#include <utility>

static inline double overlap(const std::pair<double, double> p1, const std::pair<double, double> p2) {
    if (p1.second > p2.second) {
        return p2.second - p1.first;
    } else {
        return p1.second - p2.first;
    }
}

template<typename Vector>
std::pair<double, double> project_axis(const std::vector<Vector>& poly, const Vector& axis) {
    double min = dot(axis, poly[0]);
    double max = min;

    for (int i = 1; i < poly.size(); i++) {
        double p = dot(axis, poly[i]);
        if (p < min) {
            min = p;
        } else if (p > max) {
            max = p;
        }
    }

    return std::make_pair(min, max);
}

template<typename Vector>
void append_axes(const std::vector<Vector>& poly, std::vector<Vector>& axes) {
    for (int i = 0; i < poly.size(); i++) {
        int j = (i + 1) % poly.size();
        Vector edge = poly[j] - poly[i];
        // assume counterclockwise, n is outward-pointing
        Vector n = -perp(edge);
        normalize(n);
        axes.push_back(n);
    }
}

template<typename Scalar, typename Vector>
bool detect_collision_sat(const std::vector<Vector>& poly1, const std::vector<Vector>& poly2,
                          Vector& position, Scalar& depth, Vector& normal, bool& at_p1,
                          bool& f2f, Vector& offset) {
    const double tol_alignment = depth_tolerance;

    Scalar smallest = std::numeric_limits<Scalar>::max();

    bool vert_at_p1;
    f2f = false;

    std::vector<Vector> axes1;
    append_axes(poly1, axes1);
    for (int i = 0; i < axes1.size(); i++) {
        auto p1 = project_axis(poly1, axes1[i]);
        auto p2 = project_axis(poly2, axes1[i]);

        double o = overlap(p1, p2);
        if (o < 0) {
            return false;
        } else {
            if (o < smallest) {
                normal = axes1[i];
                smallest = overlap(p1, p2);
            }
        }
    }

    std::vector<Vector> axes2;
    append_axes(poly2, axes2);
    for (int i = 0; i < axes2.size(); i++) {
        auto p1 = project_axis(poly1, axes2[i]);
        auto p2 = project_axis(poly2, axes2[i]);

        double o = overlap(p1, p2);
        if (o < 0) {
            return false;
        } else {
            if (o < smallest) {
                normal = axes2[i];
                smallest = overlap(p1, p2);
            }
        }
    }

    vert_at_p1 = false;
    auto p1 = project_axis(poly1, normal);
    auto p2 = project_axis(poly2, normal);
    if (p1.second > p2.second) {
        normal = -normal; // always i to j
        p1 = project_axis(poly1, normal);
        p2 = project_axis(poly2, normal);
    }
    int verts = 0;
    for (int i = 0; i < poly1.size(); i++) {
        double d = dot(normal, poly1[i]);
        if (d > p2.first - tol) {
            position = poly1[i];
            verts++;
        }
    }
    if (verts == 1) {
        vert_at_p1 = true;
    } else {
        for (int i = 0; i < poly2.size(); i++) {
            double d = dot(normal, poly2[i]);
            if (d < p1.second + tol) {
                position = poly2[i];
            }
        }
    }

    at_p1 = vert_at_p1;

    // take care of edge to edge case
    std::vector<Vector> candidates;
    double dist = dot(normal, position);
    for (int i = 0; i < poly1.size(); i++) {
        double d = dot(normal, poly1[i]);
        if (std::abs(d - dist) < tol_alignment) {
            candidates.push_back(poly1[i]);
        }
    }

    for (int i = 0; i < poly2.size(); i++) {
        double d = dot(normal, poly2[i]);
        if (std::abs(d - dist) < tol_alignment) {
            candidates.push_back(poly2[i]);
        }
    }

    if (candidates.size() >= 4) {
        std::sort(candidates.begin(), candidates.end(), [](const Vector& a, const Vector& b) {
            if (std::abs(a(0) - b(0)) < eps) { return a(1) > b(1); }
            else { return a(0) > b(0); }
        });
        f2f = true;
        position = (candidates[1] + candidates[2]) / 2.0; // always use the middle 2
        offset = candidates[1] - position;
    }

    depth = -smallest;

    return true;
}

#endif // SAT_H
