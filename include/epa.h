#ifndef EPA_H
#define EPA_H

#include <assert.h>

#include "vec.h"

#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Core>

// finds furthest point from poly in direction dir
template<typename Vector>
auto find_furthest(const std::vector<Vector>& poly, const Vector& dir) -> Vector {
    double max_dist = std::numeric_limits<double>::lowest();
    Vector ret = *poly.begin();
    for (auto& vert : poly) {
        double dist = dir.dot(vert);
        if (dist > max_dist) {
            ret = vert;
            max_dist = dist;
        }
    }
    return ret;
}

// finds furthest point from poly in direction dir
template<typename Vector>
auto find_furthest(const std::vector<Vector>& poly, const Vector& dir, int& id) -> Vector {
    double max_dist = std::numeric_limits<double>::lowest();
    Vector ret = *poly.begin();
    for (int i = 0; i < poly.size(); i++) {
        double dist = dir.dot(poly[i]);
        if (dist > max_dist) {
            ret = poly[i];
            max_dist = dist;
            id = i;
        }
    }
    return ret;
}

/**
 * The functions returns the triangle simplex that GJK finds if there is intersection
 * @tparam Vector
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 */
template<typename Vector>
auto find_simplex_gjk(const std::vector<Vector>& vtxsA, const std::vector<Vector>& vtxsB) -> std::vector<Vector> {
    // compute support on Minkowski difference
    auto support = [&vtxsA, &vtxsB](const Vector& dir) -> Vector {
        Vector ndir = -dir;
        return find_furthest(vtxsA, dir) - find_furthest(vtxsB, ndir);
    };

    // add first and second points of simplex
    std::vector<Vector> simplex;
    Vector init_dir = vtxsB[0] - vtxsA[0];
    simplex.push_back(support(init_dir));
    Vector O = Vector{0, 0};
    Vector d = O - simplex[0];

    Vector AB, AC, AO, ABperp, ACperp;  // all temporary variables
    while (true) {
        Vector P = support(d);
        if (P.dot(d) < 0) {
            return std::vector<Vector>();
        } else {
            simplex.push_back(P);
        }

        if (simplex.size() == 2) {  // line case
            AB = simplex[0] - simplex[1];
            AO = O - simplex[1];
            d = AO - AB.dot(AO) / AB.squaredNorm() * AB;  // ABperp
        } else if (simplex.size() == 3) {                 // triangle case
            AB = simplex[1] - simplex[2];
            AC = simplex[0] - simplex[2];
            AO = O - simplex[2];
            ABperp = -(AC - AB.dot(AC) / AB.squaredNorm() * AB);
            ACperp = -(AB - AC.dot(AB) / AC.squaredNorm() * AC);
            if (ABperp.dot(AO) > 0) {
                simplex.erase(simplex.begin());  // remove C
                d = ABperp;
            } else if (ACperp.dot(AO) > 0) {
                simplex.erase(simplex.begin() + 1);  // remove B
                d = ACperp;
            } else {
                return simplex;
            }
        }
    }
}

/**
 * @brief return normal vector towards origin of the closest edge and the starting index
 * @param simplex
 */
template<typename Vector>
auto find_closest_edge(const std::vector<Vector>& simplex) -> std::pair<int, Vector> {
    double min_dist = std::numeric_limits<double>::max();
    int min_idx;
    Vector ret_normal;

    for (int i = 0; i < simplex.size(); i++) {
        int j = (i + 1 + simplex.size()) % simplex.size();
        Vector edge = simplex[j] - simplex[i];
        Vector n{edge[1], -edge[0]};  // we know the simplex is counterclockwise, so origin is always at left
        n.normalize();
        double dist = n.dot(simplex[i]);
        if (dist < min_dist) {
            min_dist = dist;
            min_idx = i;
            ret_normal = n;
        }
    }

    return {min_idx, ret_normal};
}

/**
 * computing maximum penetration depth and its normal for the intersection of convex hulls
 * @tparam Vector Eigen::Vector2x
 * @param[out] normalA if we move all the vertices of vtxB with normalA, there is no collision
 * @param[in] vtxsA coordinates of point set A
 * @param[in] vtxsB coordinates of point set B
 * @return Direction (Vector), depth
 */
template<typename Vector>
bool penetration_epa(Vector& normalA, const std::vector<Vector>& vtxsA, const std::vector<Vector>& vtxsB, std::vector<Vector>& simplex, double tolerance) {
    // std::vector<Vector> simplex = FindSimplex_Gjk(vtxsA, vtxsB);
    if (simplex.size() == 0) {  // no intersection
        return false;
    }
    assert(simplex.size() == 3);

    // make the simplex counterclockwise
    double orientation = simplex[0](0) * (simplex[1](1) - simplex[2](1)) 
                       - simplex[0](1) * (simplex[1](0) - simplex[2](0)) 
                       + simplex[1](0) * simplex[2](1) - simplex[1](1) * simplex[2](0);
    if (orientation < 0.0) {  // clockwise
        Vector temp = simplex[2];
        simplex[2] = simplex[1];
        simplex[1] = temp;
    }

    auto support = [&vtxsA, &vtxsB](const Vector& dir) -> Vector {
        Vector ndir = -dir;
        return find_furthest(vtxsA, dir) - find_furthest(vtxsB, ndir);
    };

    while (true) {
        std::pair<int, Vector> ret = find_closest_edge(simplex);
        int v0 = ret.first;
        Vector n = ret.second;
        double dist = n.dot(simplex[v0]);

        Vector p = support(n);
        double d = p.dot(n);

        if (d - dist < tolerance) {
            normalA = n * dist;
            return true;
        } else {
            simplex.insert(simplex.begin() + v0 + 1, p);
        }
    }
}

template<typename Vector>
struct simplex_node {
    Vector p;
    int i;
    int j;
};

template<typename Vector>
Vector closest_point(const Vector& a, const Vector& b) {
    Vector n = perp(a - b);
    double dist = dot(n, a);
    Vector c = n * dist;
    if (dot(a - c, c - b) > 0) {
        return c;
    } else {
        if (norm2(a) > norm2(b)) {
            return b;
        } else {
            return a;
        }
    }
}

/**
 * @tparam Vector
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 */
template<typename Vector>
bool find_closest_gjk(const std::vector<Vector>& vtxsA, const std::vector<Vector>& vtxsB,
                      simplex_node<Vector>& ret_a, simplex_node<Vector>& ret_b, double& depth) {
    const double tol = 1E-6;

    // compute support on Minkowski difference
    auto support = [&vtxsA, &vtxsB](const Vector& dir) -> simplex_node<Vector> {
        simplex_node<Vector> ret;
        Vector ndir = -dir;
        ret.p = find_furthest(vtxsA, dir, ret.i) - find_furthest(vtxsB, ndir, ret.j);
        return ret;
    };

    // add first and second points of simplex
    std::vector<simplex_node<Vector> > simplex;
    Vector init_dir = vtxsB[0] - vtxsA[0];
    simplex.push_back(support(init_dir));
    simplex.push_back(support(-init_dir));

    Vector d = closest_point(simplex[0].p, simplex[1].p);
    while (true) {
        if (norm2(d) < tol) {
            return false;
        }
        d = -d;
        simplex_node<Vector> c = support(d);

        double dc = dot(c.p, d);
        double da = dot(simplex[0].p, d);

        if (dc - da < tol) {
            depth = norm(d);
            ret_a = simplex[0];
            ret_b = simplex[1];
            return true;
        }

        Vector p1 = closest_point(simplex[0].p, c.p);
        Vector p2 = closest_point(c.p, simplex[1].p);

        if (norm2(p1) < norm2(p2)) {
            simplex[1] = c;
            d = p1;
        } else {
            simplex[0] = c;
            d = p2;
        }
    }
}

#endif // EPA_H
