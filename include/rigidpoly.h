#ifndef RIGIDPOLY_H
#define RIGIDPOLY_H

#include "collision_parameters.h"
#include "contact.h"
#include "epa.h"
#include "sat.h"
#include "vec.h"

#include <limits>
#include <initializer_list>
#include <vector>

template<typename Scalar, typename Vector>
struct RigidPoly {
    typedef Vec3i tri;
    std::vector<Vector> verts;
    std::vector<tri> triangles;

    // rigid body related variables
    typedef Contact<Scalar, Vector> contact_type;

    Vector center;
    Scalar rot;

    // primary variables
    Scalar w; // weight
    Vector m; // momentum
    Scalar L; // angular momentum
    Scalar I; // "inertia tensor"

    Scalar w_inv;
    Scalar I_inv;

    Scalar density;

    // derived variables
    Vector v; // velocity (d)
    Scalar omega; // angular velocity (d)

    // computed variables
    Vector F; // force
    Scalar torque; // torque (perpendicular to the plane)
    std::vector<Vector> world; // verts in world coordinate

    // other parameters
    bool fixed;

    RigidPoly() 
        : verts(), fixed(false) {}
    
    RigidPoly(std::initializer_list<Vector> vertices, Scalar density)
        : verts(vertices), center({0.0, 0.0}), rot(0.0), m(0.0), L(0.0), density(density), fixed(false) {
            for (int i = 0; i < verts.size(); i++) {
                world.push_back(verts[i]);
            }
            triangulate();
            correct_center();
            init_inertia();
        }

    contact_type detect_collision(RigidPoly& other);

    void triangulate();

    void correct_center();
    void init_inertia();

    void local_to_world();
    Vector calculate_world_vel(Vector r);

    void set_fixed();
};

template<typename Scalar, typename Vector>
typename RigidPoly<Scalar, Vector>::contact_type 
RigidPoly<Scalar, Vector>::detect_collision(RigidPoly<Scalar, Vector>& other) {
    RigidPoly<Scalar, Vector>::contact_type ret;

    // GJK
    // first check for intersection
    // std::vector<Vector> simplex = find_simplex_gjk(world, other.world);
    // if (simplex.size() > 0) { // intersection case
    //     Vector normal;
    //     Scalar depth = penetration_epa(normal, world, other.world, simplex, tol);
    // } else {
    //     simplex_node<Vector> a;
    //     simplex_node<Vector> b;
    //     Scalar depth;
    //     find_closest_gjk(world, other.world, a, b, depth);        
    // }

    // SAT
    if (!detect_collision_sat(world, other.world, ret.p, ret.depth, ret.direction, ret.vert_at_i, ret.f2f, ret.f2f_offset)) {
        ret.depth = std::numeric_limits<Scalar>::max();
        return ret;
    }

    ret.v_rel = other.v - v; // translational velocity
    ret.r_i = ret.p - center;
    ret.v_rel += cross(ret.r_i, omega);
    ret.r_j = ret.p - other.center;
    ret.v_rel -= cross(ret.r_j, other.omega);

    return ret;
}

template<typename Scalar, typename Vector>
static inline Scalar triangle_area(const Vector& p0, const Vector& p1, const Vector& p2) {
    return static_cast<Scalar>(0.5) * (std::abs(p0(0) * (p1(1) - p2(1)) 
                                              + p1(0) * (p2(1) - p0(1)) 
                                              + p2(0) * (p0(1) - p1(1))));
}

template<typename Vector>
static inline Vector triangle_centroid(const Vector& p0, const Vector& p1, const Vector& p2) {
    return (p0 + p1 + p2) / 3.0;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::triangulate() {
    // using simple triangulation for now
    for (int i = 1; i < verts.size() - 1; i++) {
        triangles.push_back({0, i, i + 1});
    }
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::correct_center() {
    Scalar sum_mass = 0.0;
    Vector weighted_sum = {0.0, 0.0};
    for (auto& tri : triangles) {
        Vector tri_centroid = triangle_centroid(verts[tri(0)], verts[tri(1)], verts[tri(2)]);
        Scalar tri_mass = triangle_area<Scalar>(verts[tri(0)], verts[tri(1)], verts[tri(2)]);

        weighted_sum += tri_mass * tri_centroid;
        sum_mass += tri_mass;
    }
    Vector center_offset = weighted_sum / sum_mass;
    center_offset -= center;
    for (auto& v : verts) {
        v -= center_offset;
    }
    center += center_offset;

    w = sum_mass * density;
    w_inv = 1.0 / w;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::init_inertia() {
    Scalar mmoi = 0.0;
    // compute inertia for each triangle
    for (auto& tri : triangles) {
        Vector tri_centroid = triangle_centroid(verts[tri(0)], verts[tri(1)], verts[tri(2)]);
        Scalar tri_mass = triangle_area<Scalar>(verts[tri(0)], verts[tri(1)], verts[tri(2)]) * density;
        Vector a = verts[tri(0)] - tri_centroid;
        Vector b = verts[tri(1)] - tri_centroid;
        Vector c = verts[tri(2)] - tri_centroid;
        Scalar tri_mmoi = static_cast<Scalar>(1.0 / 6.0) * tri_mass * (dot(a, a) + 
                                                                       dot(b, b) +
                                                                       dot(c, c) +
                                                                       dot(a, b) +
                                                                       dot(b, c) +
                                                                       dot(c, a));
        mmoi += dist2(tri_centroid, center) * tri_mass + tri_mmoi;
    }

    I = mmoi;
    I_inv = 1.0 / mmoi;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::local_to_world() {
    Scalar a11 = std::cos(rot);
    Scalar a12 = -std::sin(rot);
    Scalar a21 = -a12;
    Scalar a22 = a11;

    Scalar t1;
    Scalar t2;

    for (int i = 0; i < verts.size(); i++) {
        world[i] = verts[i];
        t1 = a11 * world[i](0) + a12 * world[i](1);
        t2 = a21 * world[i](0) + a22 * world[i](1);

        world[i](0) = t1;
        world[i](1) = t2;

        world[i] += center;
    }
}

template<typename Scalar, typename Vector>
Vector RigidPoly<Scalar, Vector>::calculate_world_vel(Vector r) {
    return v - cross(r, omega);
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::set_fixed() {
    w_inv = 0;
    I_inv = 0;
    fixed = true;
}

#endif
