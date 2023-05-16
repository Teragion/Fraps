#include "util.h"
#include "util/color.h"
#include "util/heatmap.h"
#include "util/monitor.h"

#include "collision_parameters.h"
#include "poly_reader.h"
#include "rigidpoly.h"
#include "triangulate.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <limits>
#include <vector>

#include <Eigen/Dense>

#include <chrono>
#include <iostream>
#include <thread>

#include <cassert>

#include <png.h>

using contact_type = Contact<Scalar, Vector>;

Scalar dt_last = dt;

// TODO: replace this with Runge-Kutta
void ode(std::vector<RigidPoly<Scalar, Vector> >& objs, Scalar dt) {
    dt_last = dt;
    // assuming F is computed and filled
    for (auto& o : objs) {
        if (o.deleted || o.fixed) {
            continue;
        }
        if (dt > 0) {
            o.center += o.v * dt;
            o.rot += o.omega * dt;
            o.m += o.F * dt;
            o.L += o.torque * dt;
            o.v = o.m * o.w_inv;
            o.omega = o.L * o.I_inv;
        } else {
            o.m += o.F * dt;
            o.L += o.torque * dt;
            o.v = o.m * o.w_inv;
            o.omega = o.L * o.I_inv;
            o.center += o.v * dt;
            o.rot += o.omega * dt;
        }

        o.local_to_world();
    }
}

Scalar detect_contact(std::vector<RigidPoly<Scalar, Vector> >& objs, 
                      std::vector<contact_type>& impulse,
                      std::vector<contact_type>& resting,
                      Scalar& suggest) {
    for (int i = 0; i < objs.size(); i++) {
        for (int j = i + 1; j < objs.size(); j++) {
            if (objs[i].deleted || objs[j].deleted) {
                continue;
            }
            contact_type contact = objs[i].detect_collision(objs[j]);
#ifdef ADAPTIVE_STEP
            if (contact.depth < -depth_tolerance) {
                return contact.depth;
            }
#endif
            if (contact.depth < depth_tolerance) {
                Scalar v_rel_s = dot(contact.v_rel, contact.direction);
                contact.v_rel_s = v_rel_s;
#ifndef NDEBUG
                std::cout << objs[i].v << std::endl;
                std::cout << objs[i].omega << std::endl;
                std::cout << objs[j].v << std::endl;
                std::cout << objs[j].omega << std::endl;
                std::cout << contact.depth << std::endl;
                std::cout << contact.r_i << std::endl;
                std::cout << contact.v_rel << std::endl;
                std::cout << contact.v_rel_s << std::endl;
#endif
                contact.i = i;
                contact.j = j;
                if (contact.f2f) {
                    contact_type c1 = contact;
                    c1.p = c1.p + c1.f2f_offset;
                    c1.v_rel = objs[j].v - objs[i].v;
                    c1.r_i = c1.p - objs[i].center;
                    c1.v_rel -= Vector{objs[i].omega * -c1.r_i(1), objs[i].omega * c1.r_i(0)};
                    c1.r_j = c1.p - objs[j].center;
                    c1.v_rel += Vector{objs[j].omega * -c1.r_j(1), objs[j].omega * c1.r_j(0)};

                    contact_type c2 = contact;
                    c2.p = c2.p - c2.f2f_offset;
                    c2.v_rel = objs[j].v - objs[i].v;
                    c2.r_i = c2.p - objs[i].center;
                    c2.v_rel -= Vector{objs[i].omega * -c2.r_i(1), objs[i].omega * c2.r_i(0)};
                    c2.r_j = c2.p - objs[j].center;
                    c2.v_rel += Vector{objs[j].omega * -c2.r_j(1), objs[j].omega * c2.r_j(0)};

                    impulse.push_back(c1);
                    impulse.push_back(c2);
                } else {
                    impulse.push_back(contact);
                }
            }
        }
    }
    return 0;
}

void calculate_v_rel(std::vector<RigidPoly<Scalar, Vector> >& objs,
                     contact_type& contact) {
    contact.v_rel = objs[contact.j].calculate_world_vel(contact.r_j)
                    - objs[contact.i].calculate_world_vel(contact.r_i);
    contact.v_rel_s = dot(contact.v_rel, contact.direction);
}

void process_impulse(std::vector<RigidPoly<Scalar, Vector> >& objs, 
                     std::vector<contact_type>& impulse) {
#ifndef NDEBUG
    int itr = 0;
    int processed = 0;
#endif
    while (true) {
#ifndef NDEBUG
        itr++;
        std::cout << "process_impulse, itr = " << itr << std::endl;
        std::cout << "\t\tprocessed = " << processed << std::endl;
#endif
        bool change = false;
        for (auto& contact : impulse) {
            if (contact.processed) {
                continue;
            }
            calculate_v_rel(objs, contact);
            if (contact.v_rel_s > -tol_rel_v) {
                continue;
            }
            change = true;
#ifndef NDEBUG
            processed++;
#endif
            // compute impulse exchange
            Scalar numerator = -(1 + eps) * contact.v_rel_s;
            // Scalar ti = contact.r_i(0) * contact.direction(1) - contact.r_i(1) * contact.direction(0);
            // Scalar tj = contact.r_j(0) * contact.direction(1) - contact.r_j(1) * contact.direction(0);
            // Scalar denominator = 1 * objs[contact.i].w_inv + 1 * objs[contact.j].w_inv
            //                      + dot(1 * objs[contact.i].I_inv * Vector{ti * -contact.r_i(0), ti * contact.r_i(1)}
            //                          + 1 * objs[contact.j].I_inv * Vector{tj * -contact.r_j(0), tj * contact.r_j(1)},
            //                            contact.direction);
            Scalar denominator =
                objs[contact.i].w_inv + objs[contact.j].w_inv + sqr(cross(contact.r_i, contact.direction)) * objs[contact.i].I_inv + sqr(cross(contact.r_j, contact.direction)) * objs[contact.j].I_inv;
            Scalar j = numerator / denominator;
            Vector J = j * contact.direction;

//            objs[contact.i].m -= J;
//            objs[contact.i].L -= cross(contact.r_i, J);
//            if (!objs[contact.j].fixed) {
//                objs[contact.j].m += J;
//                objs[contact.j].L += cross(contact.r_j, J);
//            }

            objs[contact.i].apply_impulse(contact.r_i, -J, contact.i_s, contact.i_r);
            objs[contact.j].apply_impulse(contact.r_j, J, contact.j_s, contact.j_r);

            contact.processed = true;
        }
        if (!change) {
            break;
        }
    }
}

void find_resting(std::vector<RigidPoly<Scalar, Vector> >& objs,
                  std::vector<contact_type>& impulse,
                  std::vector<contact_type>& resting) {
    for (auto& contact : impulse) {
        contact.v_rel = objs[contact.j].calculate_world_vel(contact.r_j)
                        - objs[contact.i].calculate_world_vel(contact.r_i);
        contact.v_rel_s = dot(contact.v_rel, contact.direction);
        int i = contact.i;
        int j = contact.j;
        if (contact.v_rel_s < tol_rel_v) {
            if (contact.v_rel_s > tol_rel_v) { // moving away
                continue;
            } else if (contact.v_rel_s > -tol_rel_v) {  // resting contact
                contact.i = i;
                contact.j = j;
                if (contact.f2f) {
                    contact_type c1 = contact;
                    c1.p = c1.p + c1.f2f_offset;
                    c1.v_rel = objs[j].v - objs[i].v;
                    c1.r_i = c1.p - objs[i].center;
                    c1.v_rel -= Vector{objs[i].omega * -c1.r_i(1), objs[i].omega * c1.r_i(0)};
                    c1.r_j = c1.p - objs[j].center;
                    c1.v_rel += Vector{objs[j].omega * -c1.r_j(1), objs[j].omega * c1.r_j(0)};

                    contact_type c2 = contact;
                    c2.p = c2.p - c2.f2f_offset;
                    c2.v_rel = objs[j].v - objs[i].v;
                    c2.r_i = c2.p - objs[i].center;
                    c2.v_rel -= Vector{objs[i].omega * -c2.r_i(1), objs[i].omega * c2.r_i(0)};
                    c2.r_j = c2.p - objs[j].center;
                    c2.v_rel += Vector{objs[j].omega * -c2.r_j(1), objs[j].omega * c2.r_j(0)};

                    resting.push_back(c1);
                    resting.push_back(c2);
                } else {
                    resting.push_back(contact);
                }
            } else {
                //                assert(false);
            }
        }
    }
}


void clear_variables(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        clear_vec(o.v);
        o.omega = 0.0;

        clear_vec(o.F);
        o.torque = 0.0;
    }
}

void update_derived(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        o.v = o.m * o.w_inv;
        o.omega = o.L * o.I_inv;
    }
}

void update_last(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        o.v_last = o.v;
        o.omega_last = o.omega;
    }
}

void apply_gravity(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        if (o.fixed) { continue; }
        o.F += g * o.w;
    }
}

Eigen::VectorX<Scalar> fdirection(const Eigen::MatrixX<Scalar>& A, const Eigen::VectorX<int>& C, int d, int n) {
    Eigen::VectorX<Scalar> delta_f = Eigen::VectorX<Scalar>::Zero(n);
    delta_f(d) = 1;

    // construct matrix Acc
    int row;
    int c = C.sum();
    Eigen::MatrixX<Scalar> Acc(c, c);
    row = 0;
    for (int i = 0; i < n; i++) {
        if (C(i) == 0) {
            continue;
        }
        int col = 0;
        for (int j = 0; j < n; j++) {
            if (C(j) == 0) {
                continue;
            }
            Acc(row, col) = A(i, j);
            col++;
        }
        row++;
    }

    // construct Acd
    Eigen::VectorX<Scalar> Acd(c);
    row = 0;
    for (int i = 0; i < n; i++) {
        if (C(i) == 0) {
            continue;
        }
        Acd(row) = -A(i, d);
        row++;
    }

    // solve
    Eigen::VectorX<Scalar> x = Acc.lu().solve(Acd);

    // transfer x to delta_f
    row = 0;
    for (int i = 0; i < n; i++) {
        if (C(i) == 1) {
            delta_f(i) = x(row);
            row++;
        }
    }

    return delta_f;
}

void max_step(const Eigen::VectorX<Scalar>& f, const Eigen::VectorX<Scalar>& a, 
              const Eigen::VectorX<Scalar>& delta_f, const Eigen::VectorX<Scalar>& delta_a, 
              const Eigen::VectorX<int>& C, const Eigen::VectorX<int>& NC,
              int d, int n, Scalar& s, int& j) {
    s = std::numeric_limits<Scalar>::max();
    j = -1;
    if (delta_a(d) > 0) {
        j = d;
        s = -a(d) / delta_a(d);
    }

    // for i in C...
    Scalar s_p;
    for (int i = 0; i < n; i++) {
        if (C(i) == 0) { continue; }
        if (delta_f(i) < -tol) {
            s_p = -f(i) / delta_f(i);
            if (s_p < s) {
                s = s_p;
                j = i;
            }
        }
    }

    // for i in NC...
    for (int i = 0; i < n; i++) {
        if (NC(i) == 0) { continue; }
        if (delta_a(i) < -tol) {
            s_p = -a(i) / delta_a(i);
            if (s_p < s) {
                s = s_p;
                j = i;
            }
        }
    }
    // just happen to realize these loops can be merged, anyway...
}

void compute_forces(const Eigen::MatrixX<Scalar>& A, const Eigen::VectorX<Scalar>& b, Eigen::VectorX<Scalar>& f, int n) {
    constexpr Scalar tol = 1E-9;

    // assuming f = 0 already
    Eigen::VectorX<Scalar> a = b;

    // initialize C and NC to be empty sets
    Eigen::VectorX<int> C = Eigen::VectorX<int>::Zero(n);
    Eigen::VectorX<int> NC = Eigen::VectorX<int>::Zero(n);

    int itr = 0;

    while(true) {
        itr++;
#ifndef NDEBUG
        std::cout << "itr = " << itr << std::endl;
#endif
        // find d to drive to zero
        int d = -1;
        for (int i = 0; i < n; i++) {
            if (a(i) < -tol) {
                d = i;
                break;
            }
        }
        if (d == -1) {
            break;
        }
#ifndef NDEBUG
        std::cout << "a(d) = " << a(d) << std::endl;
#endif

        // drive to zero (d)

        // fdirection
        Eigen::VectorX<Scalar> delta_f = fdirection(A, C, d, n);
        Eigen::VectorX<Scalar> delta_a = A * delta_f;
#ifndef NDEBUG
        std::cout << "delta_f = " << delta_f << std::endl;
        std::cout << "delta_a = " << delta_a << std::endl;
#endif

        // maxstep
        Scalar s;
        int j;
        max_step(f, a, delta_f, delta_a, C, NC, d, n, s, j);

        f = f + s * delta_f;
        a = a + s * delta_a;

        if (C(j) == 1) {
            C(j) = 0;
            NC(j) = 1;
        } else if (NC(j) == 1) {
            C(j) = 1;
            NC(j) = 0;
        } else {
            C(j) = 1;
        }
    }
}

void solve_contact_forces(std::vector<RigidPoly<Scalar, Vector> >& objs, std::vector<contact_type>& resting) {
    int n_contacts = resting.size();

    if (n_contacts == 0) { return; }

#ifndef NDEBUG
    std::cout << "resting = " << n_contacts << std::endl;
#endif

    Eigen::MatrixX<Scalar> A = Eigen::MatrixX<Scalar>::Zero(n_contacts, n_contacts);
    Eigen::VectorX<Scalar> b = Eigen::VectorX<Scalar>::Zero(n_contacts);
    Eigen::VectorX<Scalar> f = Eigen::VectorX<Scalar>::Zero(n_contacts);

    // computing matrix A
    for (int i = 0; i < n_contacts; i++) {
        for (int j = 0; j < n_contacts; j++) {
            // calculating a_ij : d_i_pp's dependence on f_j
            int sign = 0;
            Scalar mass_inv; // of object of interest
            Scalar I_inv; // of object of interest
            Vector ri; // of object of interest
            Vector rj; // of object of interest

            int body1 = resting[i].i;
            int body2 = resting[i].j;

            A(i, j) = 0;

            if (resting[j].i == body1) {
                // -fjnj affecting body 1
                sign = 1;
                mass_inv = objs[body1].w_inv;
                I_inv = objs[body1].I_inv;
                ri = resting[i].r_i;
                rj = resting[j].r_i;
                if (objs[body1].fixed) { sign = 0; }
            } else if (resting[j].j == body1) {
                // fjnj affecting body 1
                sign = -1;
                mass_inv = objs[body1].w_inv;
                I_inv = objs[body1].I_inv;
                ri = resting[i].r_i;
                rj = resting[j].r_j;
                if (objs[body1].fixed) { sign = 0; }
            } if (sign != 0) {
                A(i, j) += sign * dot(resting[i].direction, resting[j].direction * mass_inv
                                                                - cross(ri, cross(rj, resting[j].direction)) * I_inv);
                sign = 0;
            }
            if (resting[j].i == body2) {
                // -fjnj affecting body 2
                sign = -1;
                mass_inv = objs[body2].w_inv;
                I_inv = objs[body2].I_inv;
                ri = resting[i].r_j;
                rj = resting[j].r_i;
                if (objs[body2].fixed) { sign = 0; }
            } else if (resting[j].j == body2) {
                // fjnj affecting body 2
                sign = 1;
                mass_inv = objs[body2].w_inv;
                I_inv = objs[body2].I_inv;
                ri = resting[i].r_j;
                rj = resting[j].r_j;
                if (objs[body2].fixed) { sign = 0; }
            } if (sign != 0) {
                A(i, j) += sign * dot(resting[i].direction, resting[j].direction * mass_inv
                                                                - cross(ri, cross(rj, resting[j].direction)) * I_inv);
            }
        }
    }

    // computing vector b
    for (int i = 0; i < n_contacts; i++) {
        int body1 = resting[i].i;
        int body2 = resting[i].j;

        // int sign;
        Scalar omega;
        if (resting[i].vert_at_i) {
            // sign = -1;
            omega = objs[body2].omega;
        } else {
            // sign = -1;
            omega = objs[body1].omega;
        }

        b(i) = -2 * dot(-cross(resting[i].direction, omega), objs[body1].v - cross(resting[i].r_i, objs[body1].omega) 
                                                                 - (objs[body2].v - cross(resting[i].r_j, objs[body2].omega)));
        Vector a1 = objs[body1].F * objs[body1].w_inv;
        Vector a2 = objs[body2].F * objs[body2].w_inv;
        b(i) += -dot(resting[i].direction, (a1 + cross(cross(resting[i].r_i, objs[body1].omega), objs[body1].omega))
                                               - (a2 + cross(cross(resting[i].r_j, objs[body2].omega), objs[body2].omega)));

        b(i) += resting[i].v_rel_s + resting[i].depth;  // penentration correction
    }

#ifndef NDEBUG
    std::cout << "finished computing A and b" << std::endl;
    std::cout << A << std::endl;
    std::cout << b << std::endl;
#endif
    compute_forces(A, b, f, n_contacts);

#ifndef NDEBUG
    std::cout << "finished computing forces" << std::endl;
    std::cout << f << std::endl;
#endif
    // apply resting contact forces
    for (int i = 0; i < n_contacts; i++) {
        int body1 = resting[i].i;
        int body2 = resting[i].j;

        // objs[body1].F -= f(i) * resting[i].direction;
        // objs[body1].torque -= cross(resting[i].r_i, f(i) * resting[i].direction);
        // objs[body2].F += f(i) * resting[i].direction;
        // objs[body2].torque += cross(resting[i].r_j, f(i) * resting[i].direction);

        objs[body1].m -= f(i) * resting[i].direction * dt;
        objs[body1].L -= cross(resting[i].r_i, f(i) * resting[i].direction) * dt;
        objs[body2].m += f(i) * resting[i].direction * dt;
        objs[body2].L += cross(resting[i].r_j, f(i) * resting[i].direction) * dt;
    }
}

bool check_fracture(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    bool is_fracture = false;
    for (auto& o : objs) {
        if (o.deleted) {
            continue;
        }
//        Scalar max_stress = 0.0;
//        for (int i = 0; i < o.sections.size(); i++) {
//            Section& s = o.sections[i];
//
//            Vector r = o.causing_impact.pos;
//
//            bool impulse_side = (dot(r, s.normal) > dot(s.p, s.normal));
//            if (dot(s.c0, s.normal) > dot(s.p, s.normal)) {
//                impulse_side = !impulse_side;
//            }
//
//            Vector a0 = o.calculate_world_vel(s.c0) - o.calculate_world_vel_last(s.c0);
//            Vector a1 = o.calculate_world_vel(s.c1) - o.calculate_world_vel_last(s.c1);
//
//            Vector f0 = a0 / o.fracture_dt / s.w0_inv;
//            Vector f1 = a1 / o.fracture_dt / s.w1_inv;
//
//            Scalar t0 = std::abs(cross(s.c0, f0));
//            Scalar t1 = std::abs(cross(s.c1, f1));
//
//            Scalar s0 = t0 * (s.length / 2) / moment_of_inertia(s.length);
//            Scalar s1 = t1 * (s.length / 2) / moment_of_inertia(s.length);
//
//            if ((s0 > o.stress_toleration && impulse_side) || (s1 > o.stress_toleration && (!impulse_side))) {
//                is_fracture = true;
//                o.to_fracture = true;
//                Scalar stress = max(s0, s1);
//                if (stress > max_stress) {
//                    max_stress = stress;
//                    o.fracture_pos = i;
//                }
//            }
//        }
        o.check_fracture();
        if (o.to_fracture) {
            is_fracture = true;

            Section& s = o.sections[o.fracture_pos];

            o.new_m0 = o.calculate_world_vel_last(s.c0) / s.w0_inv;
            o.new_m1 = o.calculate_world_vel_last(s.c1) / s.w1_inv;
            o.new_L0 = o.omega / s.I0_inv;
            o.new_L1 = o.omega / s.I1_inv;

            o.distribute_impulse();
        }
    }

    return is_fracture;
}

void rewind(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        o.v = o.v_last;
        o.omega = o.omega_last;
    }
}

/**
 * detach all bodies that is marked "detached" along "fracture_pos"
 * @param objs
 */
void detach(std::vector<RigidPoly<Scalar, Vector> >& objs) {
    unsigned int size = objs.size();
    for (int k = 0; k < size; k++) {
        if (objs[k].detached) {
            objs[k].deleted = true;
            unsigned int idx = objs.size();
            objs.emplace_back();
            objs.emplace_back();

            objs[idx].density = objs[k].density;
            objs[idx + 1].density = objs[k].density;

            auto& s0 = objs[idx];
            auto& s1 = objs[idx + 1];
            auto& section = objs[k].sections[objs[k].fracture_pos];

            int cur_side = 0;

            for (int i = 0; i < objs[k].verts.size(); i++) {
                if (cur_side == 0) {
                    s0.verts.push_back(objs[k].verts[i]);
                } else {
                    s1.verts.push_back(objs[k].verts[i]);
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

//            s0.triangulate();
//            s0.correct_center();
            s0.init_inertia();

//            s1.triangulate();
//            s1.correct_center();
            s1.init_inertia();

            if (objs[k].decomposed) {
                s0.convex_decompose();
                s1.convex_decompose();
            }

            for (int i = 0; i < s0.verts.size(); i++) {
                s0.world.push_back(s0.verts[i]);
            }

            for (int i = 0; i < s1.verts.size(); i++) {
                s1.world.push_back(s1.verts[i]);
            }

            s0.m = objs[k].new_m0;
            s1.m = objs[k].new_m1;

            s0.L = objs[k].new_L0;
            s1.L = objs[k].new_L1;

            s0.center = objs[k].center + objs[k].local_to_world(section.c0);
            s0.rot = objs[k].rot;
            s1.center = objs[k].center + objs[k].local_to_world(section.c1);
            s1.rot = objs[k].rot;

            s0.local_to_world();
            s1.local_to_world();

            // s0.enable_fracture = false;
            // s1.enable_fracture = false;
        }
    }

    auto it = std::remove_if(objs.begin(), objs.end(), [](RigidPoly<Scalar, Vector> o) { return o.deleted; });
    objs.erase(it, objs.end());
}

int main() {
    // initialize objects
    std::vector<RigidPoly<Scalar, Vector> > objs;

//    objs.push_back(RigidPoly<Scalar, Vector>({{-5.0, -1.0}, {-5.0, 1.0}, {-6.0, 1.0}, {-6.0, -1.0}}, 1.0));
//    objs.push_back(RigidPoly<Scalar, Vector>({{-5.0, 0.0}, {-5.0, 0.7}, {-5.5, 0.8}, {-5.5, 1.2}, {-5.0, 1.3},
//                                                 {-5.0, 2.0}, {-6.5, 2.0}, {-6.5, 1.3}, {-6.0, 1.2}, {-6.0, 0.8},
//                                                 {-6.5, 0.7}, {-6.5, 0.0}}, 1.0));
    // objs.push_back(RigidPoly<Scalar, Vector>({{-4.0, -1.0}, {-2, 1.0}, {-2, 3.8}, {-4, 5.8}, {-6.8, 5.8},
    //                                              {-8.8, 3.8}, {-8.8, 1.0}, {-6.8, -1.0}, {-5.8, -0.8}, {-6.6, -0.8},
    //                                              {-8.6, 1.0}, {-8.6, 3.8}, {-6.6, 5.6}, {-4.2, 5.6},
    //                                              {-2.2, 3.8}, {-2.2, 1.0}, {-4.1, -0.888}}, 1.0));

    objs.push_back(PolyReader::readSolid("/Users/teragion/Models/out.poly"));
    objs[0].enable_fracture = true;
    // objs[0].convex_decompose();
    Triangulate::triangulateSolid(objs[0]);
    objs[0].stress_toleration = 10;

//    for (int i = 0; i < 6; i++) {
//        objs.push_back(RigidPoly<Scalar, Vector>({{i - 7.5, 0.0}, {i - 7.5, 4.0}, {i - 8.0, 4.0}, {i - 8.0, 0.0}}, 4.0));
//        objs[i].enable_fracture = true;
//        objs[i].stress_toleration = 18000;
//    }

    objs.push_back(RigidPoly<Scalar, Vector>({{4.5, 0.6}, {4.5, 0.8}, {4.2, 0.7}}, 2.0));
    objs[1].enable_fracture = false;
    objs[1].m = {-5, 0};

//    objs.push_back(RigidPoly<Scalar, Vector>({{-8.5, 4.0}, {-2.0, 4.0}, {-5.25, 4.5}}, 1.0));
//    objs[7].enable_fracture = false;

//    objs.push_back(RigidPoly<Scalar, Vector>({{-5.0, -4.0}, {-5.0, -2.0}, {-6.5, -2.0}, {-6.5, -4.0}}, 1.0));
//    objs[2].enable_fracture = true;
//    objs[2].stress_toleration = 6000;
//
//    objs.push_back(RigidPoly<Scalar, Vector>({{4.5, -2.2}, {6.0, -3.0}, {6.0, -1.4}}, 1.0));
//    objs[3].enable_fracture = false;
//    objs[3].m = {-15, 0};

//    objs.push_back(RigidPoly<Scalar, Vector>({{-10.5, 0.0}, {-10.5, -0.5}, {10.5, -0.5}, {10.5, 0.0}}, 1.0));
//    objs[8].set_fixed();

//    objs.push_back(RigidPoly<Scalar, Vector>({{-7.5, 0.3}, {-7.5, 0.0}, {-6.5, 0.0}, {-6.5, 0.3}}, 1.0));
//    objs[2].set_fixed();

    // objs.push_back(RigidPoly<Scalar, Vector>({{-10, -5.0}, {-9.5, -5.0}, {-9.5, 11.0}, {-10, 11.0}}, 1.0));
    // objs[2].set_fixed();

    for (auto& obj : objs) {
        if (obj.enable_fracture) {
            obj.init_sections();
        }
    }

    Monitor mon(1024, 1024);

    Heatmap color_map;
    
    glMatrixMode(GL_PROJECTION);
    glScalef(1.0, 1.0, 1.0);

    uint8_t *pixels = new uint8_t[1024 * 1024 * 3];

    unsigned int frame_count = 0;

    const double t0 = glfwGetTime();
    double lframe = t0;
    double t = 0.0;
    Scalar suggest;
    
    bool stop_flag = false;

    while (!mon.shouldClose) {
        const double t1 = glfwGetTime();
        if (!stop_flag) {
    //        std::cout << "t: " << t << std::endl;
    //        std::cout << "t1: " << t1 << std::endl;
            if (t > t1) {
                std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(std::floor((t - t1) * 1000))));
            }

            clear_variables(objs);
            update_derived(objs);
            update_last(objs);
            apply_gravity(objs);

            std::vector<contact_type> resting;
            std::vector<contact_type> impulse;
#ifndef NDEBUG
            std::cout << "detect_contact()"<< std::endl;
#endif
            detect_contact(objs, impulse, resting, suggest);

#ifndef NDEBUG
            std::cout << "process_impulse()"<< std::endl;
#endif
            process_impulse(objs, impulse);
#ifndef NDEBUG
            std::cout << "find_resting()"<< std::endl;
#endif
            update_derived(objs);
            find_resting(objs, impulse, resting);

#ifndef NDEBUG
            std::cout << "solve_contact_forces()"<< std::endl;
#endif
            solve_contact_forces(objs, resting);
            update_derived(objs);

            if (check_fracture(objs)) {
                stop_flag = true;
                continue;
                rewind(objs);
                detach(objs);
                continue;
            }

            ode(objs, dt);
        }
        t += dt;
#ifdef ADAPTIVE_STEP
        // fine tune the step size
        bool flag = false;
        int attempt = 0;
        Scalar cur_step = dt / 2.0;
        while (true) {
            attempt++;
            if (attempt > 20) {
                return -1;
            }
            Scalar ret = detect_contact(objs, impulse, resting, suggest);
#ifndef NDEBUG
            std::cout << "depth: " << ret << std::endl;
            std::cout << "t: " << t << std::endl;
#endif
            if (ret < -depth_tolerance) {
                ode(objs, -cur_step);
                t -= cur_step;
                cur_step /= 2.0;
                flag = true;
            } else if (!flag) {
                break;
            } else if (ret > depth_tolerance) {
                ode(objs, cur_step);
                t += cur_step;
                cur_step /= 2.0;
            } else {
                break;
            }
        }
#endif

        if (t1 - lframe > draw_step) {
            frame_count++;
//            std::cout << "Frame: " << frame_count << std::endl;
            lframe = t1;
            // drawing procedure
            mon.clear();

            for (int i = 0; i < objs.size(); i++) {
                auto& c = objs[i];
                Color color;
                color_map.getColorAtValue(static_cast<float>(i) / static_cast<float>(objs.size()), color);
                if (objs[i].triangulated) {
                    // for (int j = 0; j < objs[i].convex_polys_world.size(); j++) {
                    //     mon.drawPoly(objs[i].convex_polys_world[j], 1, white);
                    // }

                    for (int j = 0; j < objs[i].triangles.size(); j++) {
                        const auto& tri = objs[i].triangles[j];
                        mon.drawTri(objs[i].inner_world[tri(0)], 
                                    objs[i].inner_world[tri(1)],
                                    objs[i].inner_world[tri(2)], white);
                    }
                } else {
                    mon.drawPoly(objs[i].world, 1, color);
                }
            }

            Scalar max_r;
            for (int i = 0; i < objs[0].sections.size(); i++) {
                if ((objs[0].sections[i].ratio) > max_r) {
                    max_r = objs[0].sections[i].ratio;
                }
            }

            for (int i = 0; i < objs[0].sections.size(); i++) {
                Color color;
                color_map.getColorAtValue(static_cast<float>(objs[0].sections[i].ratio) / max_r, color);
                mon.drawLine(objs[0].center + objs[0].sections[i].j0.p, objs[0].center + objs[0].sections[i].j1.p, 1, color);
            }

            // write to png
//            glReadPixels(0, 0, 1024, 1024, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);
//            for (int j = 0; j * 2 < 1024; ++j) {
//                int x = j * 1024 * 3;
//                int y = (1024 - 1 - j) * 1024 * 3;
//                for (int i = 1024 * 3; i > 0; --i) {
//                    uint8_t tmp = pixels[x];
//                    pixels[x] = pixels[y];
//                    pixels[y] = tmp;
//                    ++x;
//                    ++y;
//                }
//            }
//            png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
//            png_infop info = png_create_info_struct(png);
//            if (!info) {
//                std::cout << "unable to create write struct" << std::endl;
//                png_destroy_write_struct(&png, &info);
//                return -1;
//            }
//            std::string s = "output/Frame" + std::to_string(frame_count) + ".png";
//            FILE *fp = fopen(s.c_str(), "wb");
//            if (!fp) {
//                std::cout << "unable to open output folder" << std::endl;
//                png_destroy_write_struct(&png, &info);
//                return -1;
//            }
//
//            png_init_io(png, fp);
//            png_set_IHDR(png, info, 1024, 1024, 8 , PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
//                         PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
//            png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
//            if (!palette) {
//                fclose(fp);
//                png_destroy_write_struct(&png, &info);
//                return -1;
//            }
//            png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
//            png_write_info(png, info);
//            png_set_packing(png);
//
//            png_bytepp rows = (png_bytepp)png_malloc(png, 1024 * sizeof(png_bytep));
//            for (int i = 0; i < 1024; ++i)
//                rows[i] = (png_bytep)(pixels + (1024 - i - 1) * 1024 * 3);
//
//            png_write_image(png, rows);
//            png_write_end(png, info);
//            png_free(png, palette);
//            png_destroy_write_struct(&png, &info);
//
//            fclose(fp);
//            delete[] rows;

            mon.refresh();
        }
    }
}
