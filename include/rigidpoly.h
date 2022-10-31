#ifndef RIGIDPOLY_H
#define RIGIDPOLY_H

#include "collision_parameters.h"
#include "contact.h"
#include "epa.h"
#include "force.h"
#include "polygon.h"
#include "sat.h"
#include "section.h"
#include "vec.h"

#include <limits>
#include <list>
#include <initializer_list>
#include <vector>

template<typename Scalar, typename Vector>
struct RigidPoly {
    typedef Vec3i tri;
    std::vector<Vector> verts;
    std::vector<tri> triangles;

    bool decomposed = false;
    typedef std::vector<Vector> Shape;
    std::vector<Shape> convex_polys;
    std::vector<Shape> convex_polys_world;
    std::vector<Vector> steinerPoints, reflexVertices; // auxillary arrays

    // rigid body related variables
    typedef Contact<Scalar, Vector> contact_type;

    Vector center;
    Scalar rot;

    Vector center_last;
    Scalar rot_last;

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

    Vector v_last;
    Scalar omega_last;

    // computed variables
    Vector F; // force
    Scalar torque; // torque (perpendicular to the plane)
    std::vector<Vector> world; // verts in world coordinate

    // fracture related varibles
    std::vector<Section> sections;
    Scalar fracture_dt = 0.01;
    Scalar stress_toleration = 2.0;
    bool to_fracture = false;
    int fracture_pos = -1;
    Impact<Scalar, Vector> causing_impact;
    bool first_impact = true;
    bool detached = false;
    Vector new_m0;
    Vector new_m1;
    Scalar new_L0;
    Scalar new_L1;
    std::list<Impact<Scalar, Vector> > impacts;

    // other parameters
    bool fixed = false;
    bool enable_fracture = true;
    bool deleted = false;

    RigidPoly() 
        : verts(), fixed(false) {}
    
    RigidPoly(std::initializer_list<Vector> vertices, Scalar density, bool init_world = true)
        : verts(vertices), center({0.0, 0.0}), rot(0.0), m(0.0), L(0.0), density(density), fixed(false), enable_fracture(false) {
        if (init_world) {
            for (int i = 0; i < verts.size(); i++) {
                world.push_back(verts[i]);
            }
        }
//        triangulate();
//        correct_center();
        init_inertia();
    }

    contact_type detect_collision(RigidPoly& other);

    void triangulate();
    void decompose(std::vector<Vector>& poly);
    void convex_decompose();

    void correct_center();
    void init_inertia();

    void calculate_section_mass_mmoi(Section& section);
    void init_sections();

    void apply_impulse(Vector r, Vector j, int i_s, Scalar i_r);
    bool check_fracture();
    void distribute_impulse();
    bool find_impulse_side(Scalar i_r, int i_s, int id);

    void local_to_world();
    Vector calculate_world_vel(Vector r);
    Vector calculate_world_vel_last(Vector r);

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
    if (!decomposed) {
        if (!other.decomposed) {
            if (!detect_collision_sat(world, other.world, ret.p, ret.depth, ret.direction, ret.vert_at_i, ret.f2f, ret.f2f_offset)) {
                ret.depth = std::numeric_limits<Scalar>::max();
                return ret;
            }
        } else {
            bool detected = false;
            for (int i = 0; i < other.convex_polys.size(); i++) {
                if (detect_collision_sat(world, other.convex_polys_world[i], ret.p, ret.depth, ret.direction, ret.vert_at_i, ret.f2f, ret.f2f_offset)) {
                    detected = true;
                    break;
                }
            }
            if (!detected) {
                ret.depth = std::numeric_limits<Scalar>::max();
                return ret;
            }
        }
    } else {
        bool detected = false;
        for (int j = 0; j < convex_polys_world.size(); j++) {
            if (!other.decomposed) {
                if (detect_collision_sat(convex_polys_world[j], other.world, ret.p, ret.depth, ret.direction, ret.vert_at_i, ret.f2f, ret.f2f_offset)) {
                    detected = true;
                }
            } else {
                for (int i = 0; i < other.convex_polys.size(); i++) {
                    if (detect_collision_sat(convex_polys_world[j], other.convex_polys_world[i], ret.p, ret.depth, ret.direction, ret.vert_at_i, ret.f2f, ret.f2f_offset)) {
                        detected = true;
                        break;
                    }
                }
            }
            if (detected) {
                break;
            }
        }
        if (!detected) {
            ret.depth = std::numeric_limits<Scalar>::max();
            return ret;
        }
    }

    ret.v_rel = other.v - v; // translational velocity
    ret.r_i = ret.p - center;
    ret.v_rel += cross(ret.r_i, omega);
    ret.r_j = ret.p - other.center;
    ret.v_rel -= cross(ret.r_j, other.omega);

    // find indices
    for (int i = 0; i < world.size(); i++) {
        Vector a = world[i] - center;
        Vector b = at(world, i + 1) - center;
        if (collinear(a, b, ret.r_i)) {
            ret.i_s = i;
            ret.i_e = (i + 1) % world.size();

            ret.i_r = norm(ret.r_i - a) / norm(b - a);
            if (ret.i_r > 1) {
                continue;
            }
            break;
        }
    }

    for (int i = 0; i < other.world.size(); i++) {
        Vector a = other.world[i] - center;
        Vector b = at(other.world, i + 1) - other.center;
        if (collinear(a, b, ret.r_j)) {
            ret.j_s = i;
            ret.j_e = (i + 1) % other.world.size();

            ret.j_r = norm(ret.r_j - a) / norm(b - a);
            if (ret.j_r > 1) {
                continue;
            }
            break;
        }
    }

    return ret;
}

template<typename Scalar, typename Vector>
static inline Scalar triangle_area(const Vector& p0, const Vector& p1, const Vector& p2) {
    return static_cast<Scalar>(0.5) * (std::abs(p0(0) * (p1(1) - p2(1)) 
                                              + p1(0) * (p2(1) - p0(1)) 
                                              + p2(0) * (p0(1) - p1(1))));
}

template<typename Scalar, typename Vector>
static inline Scalar oriented_triangle_area(const Vector& p0, const Vector& p1) {
    return static_cast<Scalar>(0.5) * (p0(0) * p1(1) - p0(1) * p1(0));
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

// helper functions for convex decomposition
inline bool is_reflex(std::vector<Vector> &poly, int i) {
    return right(at(poly, i - 1), at(poly, i), at(poly, i + 1));
}

void makeCCW(std::vector<Vector> &poly) {
    int br = 0;

    // find bottom right point
    for (int i = 1; i < poly.size(); ++i) {
        if (poly[i](1) < poly[br](1) || (poly[i](1) == poly[br](1) && poly[i](0) > poly[br](0))) {
            br = i;
        }
    }

    // reverse poly if clockwise
    if (!left(at(poly, br - 1), at(poly, br), at(poly, br + 1))) {
        std::reverse(poly.begin(), poly.end());
    }
}

inline Vector intersection(const Vector &p1, const Vector &p2, const Vector &q1, const Vector &q2) {
    Vector i;
    Scalar a1, b1, c1, a2, b2, c2, det;
    a1 = p2(1) - p1(1);
    b1 = p1(0) - p2(0);
    c1 = a1 * p1(0) + b1 * p1(1);
    a2 = q2(1) - q1(1);
    b2 = q1(0) - q2(0);
    c2 = a2 * q1(0) + b2 * q1(1);
    det = a1 * b2 - a2 * b1;
    if (std::abs(det) > alignment) { // lines are not parallel
        i(0) = (b2 * c1 - b1 * c2) / det;
        i(1) = (a1 * c2 - a2 * c1) / det;
    }
    return i;
}

/**
 * See https://mpen.ca/406/bayazit
 * @tparam Scalar
 * @tparam Vector
 * @param poly
 */
template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::decompose(std::vector<Vector>& poly) {
    Vector upperInt, lowerInt, p, closestVert;
    Scalar upperDist, lowerDist, d, closestDist;
    int upperIndex, lowerIndex, closestIndex = 0;
    std::vector<Vector> lowerPoly, upperPoly;
#ifndef NDEBUG
    printf("decompose\n");
#endif
    for (int i = 0; i < poly.size(); ++i) {
        if (is_reflex(poly, i)) {
            reflexVertices.push_back(poly[i]);
            upperDist = lowerDist = std::numeric_limits<Scalar>::max();
            for (int j = 0; j < poly.size(); ++j) {
                if (left(at(poly, i - 1), at(poly, i), at(poly, j))
                    && rightOn(at(poly, i - 1), at(poly, i), at(poly, j - 1))) { // if line intersects with an edge
                    p = intersection(at(poly, i - 1), at(poly, i), at(poly, j), at(poly, j - 1)); // find the point of intersection
                    if (right(at(poly, i + 1), at(poly, i), p)) { // make sure it's inside the poly
                        d = dist2(poly[i], p);
                        if (d < lowerDist) { // keep only the closest intersection
                            lowerDist = d;
                            lowerInt = p;
                            lowerIndex = j;
                        }
                    }
                }
                if (left(at(poly, i + 1), at(poly, i), at(poly, j + 1))
                    && rightOn(at(poly, i + 1), at(poly, i), at(poly, j))) {
                    p = intersection(at(poly, i + 1), at(poly, i), at(poly, j), at(poly, j + 1));
                    if (left(at(poly, i - 1), at(poly, i), p)) {
                        d = dist2(poly[i], p);
                        if (d < upperDist) {
                            upperDist = d;
                            upperInt = p;
                            upperIndex = j;
                        }
                    }
                }
            }

            // if there are no vertices to connect to, choose a point in the middle
            if (lowerIndex == (upperIndex + 1) % poly.size()) {
#ifndef NDEBUG
                printf("Case 1: Vertex(%d), lowerIndex(%d), upperIndex(%d), poly.size(%d)\n", i, lowerIndex, upperIndex, (int) poly.size());
#endif
                p(0) = (lowerInt(0) + upperInt(0)) / 2;
                p(1) = (lowerInt(1) + upperInt(1)) / 2;
                steinerPoints.push_back(p);

                if (i < upperIndex) {
                    lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.begin() + upperIndex + 1);
                    lowerPoly.push_back(p);
                    upperPoly.push_back(p);
                    if (lowerIndex != 0) upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.end());
                    upperPoly.insert(upperPoly.end(), poly.begin(), poly.begin() + i + 1);
                } else {
                    if (i != 0) lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.end());
                    lowerPoly.insert(lowerPoly.end(), poly.begin(), poly.begin() + upperIndex + 1);
                    lowerPoly.push_back(p);
                    upperPoly.push_back(p);
                    upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.begin() + i + 1);
                }
            } else {
                // connect to the closest point within the triangle
#ifndef NDEBUG
                printf("Case 2: Vertex(%d), closestIndex(%d), poly.size(%d)\n", i, closestIndex, (int) poly.size());
#endif

                if (lowerIndex > upperIndex) {
                    upperIndex += poly.size();
                }
                closestDist = std::numeric_limits<Scalar>::max();
                for (int j = lowerIndex; j <= upperIndex; ++j) {
                    if (leftOn(at(poly, i - 1), at(poly, i), at(poly, j))
                        && rightOn(at(poly, i + 1), at(poly, i), at(poly, j))) {
                        d = dist2(at(poly, i), at(poly, j));
                        if (d < closestDist) {
                            closestDist = d;
                            closestVert = at(poly, j);
                            closestIndex = j % poly.size();
                        }
                    }
                }

                if (i < closestIndex) {
                    lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.begin() + closestIndex + 1);
                    if (closestIndex != 0) upperPoly.insert(upperPoly.end(), poly.begin() + closestIndex, poly.end());
                    upperPoly.insert(upperPoly.end(), poly.begin(), poly.begin() + i + 1);
                } else {
                    if (i != 0) lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.end());
                    lowerPoly.insert(lowerPoly.end(), poly.begin(), poly.begin() + closestIndex + 1);
                    upperPoly.insert(upperPoly.end(), poly.begin() + closestIndex, poly.begin() + i + 1);
                }
            }

            // solve the smallest poly first
            if (lowerPoly.size() < upperPoly.size()) {
                decompose(lowerPoly);
                decompose(upperPoly);
            } else {
                decompose(upperPoly);
                decompose(lowerPoly);
            }
            return;
        }
    }
    convex_polys.push_back(poly);
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::convex_decompose() {
    decompose(verts);
    decomposed = true;
    for (int i = 0; i < convex_polys.size(); i++) {
        convex_polys_world.push_back(std::vector<Vector>());
        for (int j = 0; j < convex_polys[i].size(); j++) {
            convex_polys_world[i].push_back(convex_polys[i][j]);
        }
    }
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::correct_center() {
    Scalar sum_mass = 0.0;
    Vector weighted_sum = {0.0, 0.0};
//    for (auto& tri : triangles) {
//        Vector tri_centroid = triangle_centroid(verts[tri(0)], verts[tri(1)], verts[tri(2)]);
//        Scalar tri_mass = triangle_area<Scalar>(verts[tri(0)], verts[tri(1)], verts[tri(2)]);
//
//        weighted_sum += tri_mass * tri_centroid;
//        sum_mass += tri_mass;
//    }
    Vector O{0, 0};
    for (int i = 0 ; i < verts.size(); i++) {
        int j = (i + 1) % verts.size();
        Vector tri_centroid = triangle_centroid(O, verts[i], verts[j]);
        Scalar tri_mass = oriented_triangle_area<Scalar>(verts[i], verts[j]);

        weighted_sum += tri_mass * tri_centroid;
        sum_mass += tri_mass;
    }
    Vector center_offset = weighted_sum / sum_mass;
//    center_offset -= center;
    for (auto& v : verts) {
        v -= center_offset;
    }
    center = center_offset;

    w = sum_mass * density;
    w_inv = 1.0 / w;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::init_inertia() {
    Scalar mmoi = 0.0;
    // compute inertia for each triangle
//    for (auto& tri : triangles) {
//        Vector tri_centroid = triangle_centroid(verts[tri(0)], verts[tri(1)], verts[tri(2)]);
//        Scalar tri_mass = triangle_area<Scalar>(verts[tri(0)], verts[tri(1)], verts[tri(2)]) * density;
//        Vector a = verts[tri(0)] - tri_centroid;
//        Vector b = verts[tri(1)] - tri_centroid;
//        Vector c = verts[tri(2)] - tri_centroid;
//        Scalar tri_mmoi = static_cast<Scalar>(1.0 / 6.0) * tri_mass * (dot(a, a) +
//                                                                       dot(b, b) +
//                                                                       dot(c, c) +
//                                                                       dot(a, b) +
//                                                                       dot(b, c) +
//                                                                       dot(c, a));
//        mmoi += dist2(tri_centroid, center) * tri_mass + tri_mmoi;
//    }
    Scalar sum_mass = 0.0;
    Vector weighted_sum = {0.0, 0.0};
    Scalar sum_mmoi = 0.0;
    Vector O{0, 0};
    for (int i = 0 ; i < verts.size(); i++) {
        int j = (i + 1) % verts.size();
        Vector tri_centroid = triangle_centroid(O, verts[i], verts[j]);
        Scalar tri_mass = oriented_triangle_area<Scalar>(verts[i], verts[j]);
        Scalar tri_mmoi = tri_mass * (dot(verts[i], verts[i]) + dot(verts[j], verts[j]) + dot(verts[i], verts[j])) / 6.0;

        weighted_sum += tri_centroid * tri_mass;
        sum_mass += tri_mass;
        sum_mmoi += tri_mmoi;
    }
    Vector center_offset = weighted_sum / sum_mass;
    for (auto& v : verts) {
        v -= center_offset;
    }
    center = center_offset;

    w = sum_mass * density;
    w_inv = 1.0 / w;

    sum_mmoi *= density;
    mmoi = sum_mmoi - w * dot(center_offset, center_offset);
    I = mmoi;
    I_inv = 1.0 / I;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::calculate_section_mass_mmoi(Section& section){
    // construct 2 sub bodies
    RigidPoly<Scalar, Vector> s0;
    s0.density = density;
    RigidPoly<Scalar, Vector> s1;
    s1.density = density;

    // check arrangement
    assert((section.j0.idx_v1 > section.j0.idx_v0) || ((section.j0.idx_v0 == verts.size() - 1) && (section.j0.idx_v1 == 0)));
    assert((section.j1.idx_v1 > section.j1.idx_v0) || ((section.j1.idx_v0 == verts.size() - 1) && (section.j1.idx_v1 == 0)));

    int cur_side = 0;

    for (int i = 0; i < verts.size(); i++) {
        if (cur_side == 0) {
            s0.verts.push_back(verts[i]);
        } else {
            s1.verts.push_back(verts[i]);
        }
        if (i == section.j0.idx_v0) {
            s0.verts.push_back(section.j0.p);
            s1.verts.push_back(section.j0.p);
            cur_side = 1 - cur_side;
        } else if (i == section.j1.idx_v0) {
            s0.verts.push_back(section.j1.p);
            s1.verts.push_back(section.j1.p);
            cur_side = 1 - cur_side;
        }
    }

//    s0.triangulate();
//    s0.correct_center();
    s0.init_inertia();

//    s1.triangulate();
//    s1.correct_center();
    s1.init_inertia();

    section.c0 = s0.center;
    section.w0 = s0.w;
    section.w0_inv = s0.w_inv;
    section.I0 = s0.I;
    section.I0_inv = s0.I_inv;

    section.c1 = s1.center;
    section.w1 = s1.w;
    section.w1_inv = s1.w_inv;
    section.I1 = s1.I;
    section.I1_inv = s1.I_inv;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::init_sections() {
    sections = find_cross_sections(verts);
    for (auto& s : sections) {
        calculate_section_mass_mmoi(s);
    }
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::apply_impulse(Vector r, Vector j, int i_s, Scalar i_r) {
    if (fixed) {
        return;
    }

    m += j;
    L += cross(r, j);

    impacts.push_back({r, j, i_s, i_r});

//    if (!to_fracture) {
//        m += j;
//        L += cross(r, j);
//    } else {
//        Section& s = sections[fracture_pos];
//        bool impulse_side = (dot(r, s.normal) > dot(s.p, s.normal));
//        if (dot(s.c0, s.normal) > dot(s.p, s.normal)) {
//            impulse_side = !impulse_side;
//        }
//        // impulse_side = if impulse at c1 side
//
//        v_last = v;
//        omega_last = omega;
//
//        m += j;
//        L += cross(r, j);
//
//        v = m * w_inv;
//        omega = L * I_inv;
//
//        Vector a0 = calculate_world_vel(s.c0) - calculate_world_vel_last(s.c0);
//        Vector a1 = calculate_world_vel(s.c1) - calculate_world_vel_last(s.c1);
//
//        Vector f0 = a0 / fracture_dt * s.w0;
//        Vector f1 = a1 / fracture_dt * s.w1;
//
//        Scalar t0 = std::abs(cross(s.c0, f0));
//        Scalar t1 = std::abs(cross(s.c1, f1));
//
//        Scalar s0 = t0 * (s.length / 2) / moment_of_inertia(s.length);
//        Scalar s1 = t1 * (s.length / 2) / moment_of_inertia(s.length);
//
//        if (impulse_side && (s0 > stress_toleration)) {
//            // all remaining impulse to s1
//            new_L0 += stress_toleration / (s.length / 2) * moment_of_inertia(s.length);
//            new_m0 += (calculate_world_vel_last(s.c0) + a0 * stress_toleration / s0) * s.w0;
//
//            j -= new_m0;
//            new_m1 += j;
//            new_L1 += cross(r, j);
//
//            detached = true;
//        } else if ((!impulse_side) && (s1 > stress_toleration)) {
//            // all remaining impulse to s0
//            new_L1 += stress_toleration / (s.length / 2) * moment_of_inertia(s.length);
//            new_m1 += (calculate_world_vel_last(s.c1) + a1 * stress_toleration / s1) * s.w1;
//
//            j -= new_m1;
//            new_m0 += j;
//            new_L0 += cross(r, j);
//
//            detached = true;
//        } else {
//            // evenly distribute
//            // TODO: is this ok?
//            new_m0 += j * (s.w0 / w);
//            new_m1 += j * (s.w1 / w);
//
//            Scalar t = cross(r, j);
//            new_L0 += t * (s.w0 / w);
//            new_L1 += t * (s.w1 / w);
//        }
//
//        v = v_last;
//        omega = omega_last;
//
//        m = v_last * w;
//        L = omega_last * I;
//    }
}

template<typename Scalar, typename Vector>
bool RigidPoly<Scalar, Vector>::check_fracture() {
    Scalar max_stress = 0.0;
    fracture_pos = -1;
    to_fracture = false;
    if (impacts.size() == 0) {
        return false;
    }

    for (int idx = 0; idx < sections.size(); idx++) {
        auto& s = sections[idx];
        Scalar t0 = 0.0;
        Scalar t1 = 0.0;

        for (auto& i : impacts) {
//            bool impulse_side = (dot(i.pos, s.normal) > dot(s.p, s.normal));
//            //        if (dot(s.c0, s.normal) > dot(s.p, s.normal)) {
//            //            impulse_side = !impulse_side;
//            //        }
//            if (dist2(i.pos, s.c0) > dist2(i.pos, s.c1)) {
//                impulse_side = true;
//            } else {
//                impulse_side = false;
//            }
            bool impulse_side = find_impulse_side(i.i_r, i.i_s, idx);
            // impulse_side = if impulse at c1 side

            Scalar domega = cross(i.pos, i.J);

            if (impulse_side) {
                // calculating torque transfer to side 0
                t0 += domega * s.I0;
//                t0 += domega * s.w0 * norm2(s.c0);
            } else {
                t1 += domega * s.I1;
//                t1 += domega * s.w1 * norm2(s.c1);
            }
            std::cout << impulse_side << std::endl;
        }
        if ((t0 > t1 && s.w0 > s.w1 * 2) || (t0 < t1 && s.w1 > s.w0 * 2)) {
            continue;
        }
        // cancelling
        Scalar transfer = std::abs(t0 - t1);
        Scalar stress = (transfer * (s.length / 2) / moment_of_inertia(s.length)) / fracture_dt;
        std::cout << stress << std::endl;

        if (stress > stress_toleration) {
            to_fracture = true;
            if (stress > max_stress) {
                max_stress = stress;
                fracture_pos = idx;
            }
        }
    }

    std::cout << max_stress << std::endl;

    return to_fracture;
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::distribute_impulse() {
    assert(to_fracture);
    assert(fracture_pos >= 0);

    auto& s = sections[fracture_pos];
    Scalar t0 = 0.0;
    Scalar t1 = 0.0;

    Scalar max_transfer = (stress_toleration * moment_of_inertia(s.length) / (s.length / 2)) * fracture_dt;

    for (auto& i : impacts) {
//        bool impulse_side = (dot(i.pos, s.normal) > dot(s.p, s.normal));
//        if (dot(s.c0, s.normal) > dot(s.p, s.normal)) {
//            impulse_side = !impulse_side;
//        }
//        if (dist2(i.pos, s.c0) > dist2(i.pos, s.c1)) {
//            impulse_side = true;
//        } else {
//            impulse_side = false;
//        }
        bool impulse_side = find_impulse_side(i.i_r, i.i_s, fracture_pos);
        // impulse_side = if impulse at c1 side

        if (!detached) {
            Scalar domega = cross(i.pos, i.J);
            Scalar t;

            if (impulse_side) {
                // calculating torque transfer to side 0
                t = domega * s.I0;
//                t += domega * norm2(s.c0) * s.w0;
                t0 += t;
            } else {
                t = domega * s.I1;
//                t += domega * norm2(s.c1) * s.w1;
                t1 += t;
            }
            Scalar transfer = std::abs(t0 - t1);
            if (transfer > max_transfer) {
                detached = true;
                // first distribute
                Scalar ratio = 1 - std::abs((transfer - max_transfer) / t);
                new_m0 += i.J * (s.w0 / w) * ratio;
                new_m1 += i.J * (s.w1 / w) * ratio;

                new_L0 += domega * s.I0 * ratio;
                new_L1 += domega * s.I1 * ratio;

            } else {
                // evenly distribute
                new_m0 += i.J * (s.w0 / w);
                new_m1 += i.J * (s.w1 / w);

                new_L0 += domega * s.I0;
                new_L1 += domega * s.I1;
            }
        } else {
            if (impulse_side) {
                Vector r = i.pos - s.c1;
                new_m1 += i.J;
                new_L1 += cross(r, i.J);
            } else {
                Vector r = i.pos - s.c0;
                new_m0 += i.J;
                new_L0 += cross(r, i.J);
            }
        }
    }
}

template<typename Scalar, typename Vector>
bool RigidPoly<Scalar, Vector>::find_impulse_side(Scalar i_r, int i_s, int id) {
    auto& s = sections[id];
    int cur_side = 0;
    for (int i = 0; i < verts.size(); i++) {
        if (i == s.j0.idx_v0) {
            if (i == i_s) {
                Scalar p_r = norm(s.j0.p - verts[i]) / norm(at(verts, i + 1) - verts[i]);
                if (p_r > i_r) {
                    return cur_side == 1;
                } else {
                    return cur_side != 1;
                }
            }
            cur_side = 1 - cur_side;
        } else if (i == s.j1.idx_v0) {
            if (i == i_s) {
                Scalar p_r = norm(s.j1.p - verts[i]) / norm(at(verts, i + 1) - verts[i]);
                if (p_r > i_r) {
                    return cur_side == 1;
                } else {
                    return cur_side != 1;
                }
            }
            cur_side = 1 - cur_side;
        }
        if (cur_side == 0) {
            if (i == i_s) { return false; }
        } else {
            if (i == i_s) { return true; }
        }
    }
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

    if (decomposed) {
        for (int i = 0; i < convex_polys.size(); i++) {
            for (int j = 0; j < convex_polys[i].size(); j++) {
                convex_polys_world[i][j] = convex_polys[i][j];
                auto& p = convex_polys_world[i][j];
                t1 = a11 * p(0) + a12 * p(1);
                t2 = a21 * p(0) + a22 * p(1);

                p(0) = t1;
                p(1) = t2;

                p += center;
            }
        }
    }
}

template<typename Scalar, typename Vector>
Vector RigidPoly<Scalar, Vector>::calculate_world_vel(Vector r) {
    return v - cross(r, omega);
}

template<typename Scalar, typename Vector>
Vector RigidPoly<Scalar, Vector>::calculate_world_vel_last(Vector r) {
    return v_last - cross(r, omega_last);
}

template<typename Scalar, typename Vector>
void RigidPoly<Scalar, Vector>::set_fixed() {
    w_inv = 0;
    I_inv = 0;
    fixed = true;
}

#endif // RIGIDPOLY_H
