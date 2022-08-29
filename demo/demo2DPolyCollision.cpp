#include "util.h"
#include "util/color.h"
#include "util/heatmap.h"
#include "util/monitor.h"

#include "collision_parameters.h"
#include "rigidpoly.h"
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

using contact_type = Contact<Scalar, Vector>;

// TODO: replace this with Runge-Kutta
void ode(std::vector<RigidPoly<Scalar, Vector> >& objs, Scalar dt) {
    // assuming F is computed and filled
    for (auto& o : objs) {
        if (o.fixed) {
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
            contact_type contact = objs[i].detect_collision(objs[j]);
            if (contact.depth < -depth_tolerance) {
                return contact.depth;
            }
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
//                if (contact.f2f) {
//                    contact_type c1 = contact;
//                    c1.p = c1.p + c1.f2f_offset;
//                    c1.v_rel = objs[j].v - objs[i].v;
//                    c1.r_i = c1.p - objs[i].center;
//                    c1.v_rel -= Vector{objs[i].omega * -c1.r_i(1), objs[i].omega * c1.r_i(0)};
//                    c1.r_j = c1.p - objs[j].center;
//                    c1.v_rel += Vector{objs[j].omega * -c1.r_j(1), objs[j].omega * c1.r_j(0)};
//
//                    contact_type c2 = contact;
//                    c2.p = c2.p - c2.f2f_offset;
//                    c2.v_rel = objs[j].v - objs[i].v;
//                    c2.r_i = c2.p - objs[i].center;
//                    c2.v_rel -= Vector{objs[i].omega * -c2.r_i(1), objs[i].omega * c2.r_i(0)};
//                    c2.r_j = c2.p - objs[j].center;
//                    c2.v_rel += Vector{objs[j].omega * -c2.r_j(1), objs[j].omega * c2.r_j(0)};
//
//                    impulse.push_back(c1);
//                    impulse.push_back(c2);
//                } else {
                    impulse.push_back(contact);
//                }
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
        itr++;
#ifndef NDEBUG
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
            objs[contact.i].m -= J;
            objs[contact.i].L -= cross(contact.r_i, J);
            if (!objs[contact.j].fixed) {
                objs[contact.j].m += J;
                objs[contact.j].L += cross(contact.r_j, J);
            }
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

int main() {
    // initialize objects
    std::vector<RigidPoly<Scalar, Vector> > objs;

    objs.push_back(RigidPoly<Scalar, Vector>({{1.0, 5.5}, {1.0, 7.0}, {-1.0, 7.0}, {-1.0, 5.5}}, 1.0));
//    objs[0].set_fixed();

    objs.push_back(RigidPoly<Scalar, Vector>({{1.5, 2.5}, {1.5, 4.0}, {-0.5, 3.0}}, 1.0));
    // objs[1].set_fixed();

    objs.push_back(RigidPoly<Scalar, Vector>({{-2.5, 0.5}, {-2.5, -0.5}, {2.5, -0.5}, {2.5, 0.5}}, 1.0));
    objs[2].set_fixed();

    objs.push_back(RigidPoly<Scalar, Vector>({{2.7, 1.5}, {2.7, -0.5}, {2.9, -0.5}, {2.9, 1.5}}, 1.0));
    objs[3].set_fixed();

    objs.push_back(RigidPoly<Scalar, Vector>({{-2.9, 1.5}, {-2.9, -0.5}, {-2.7, -0.5}, {-2.7, 1.5}}, 1.0));
    objs[4].set_fixed();

    Monitor mon(1024, 1024);

    Heatmap color_map;
    
    glMatrixMode(GL_PROJECTION);
    glScalef(0.05, 0.05, 1.0);

    unsigned int frame_count = 0;

    const double t0 = glfwGetTime();
    double lframe = t0;
    double t = 0.0;
    Scalar suggest;

    while (!mon.shouldClose) {
        std::cout << "t: " << t << std::endl;
        const double t1 = glfwGetTime();
        std::cout << "t1: " << t1 << std::endl;
        if (t > t1) {
            std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(std::floor((t - t1) * 1000))));
        }

        clear_variables(objs);
        apply_gravity(objs);
        update_derived(objs);

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

        // fine tune the step size
        Scalar cur_step = dt / 2.0;
        ode(objs, dt);
        t += dt;
        bool flag = false;
        int attempt = 0;
        while (true) {
            attempt++;
            if (attempt > 20) {
                return -1;
            }
            Scalar ret = detect_contact(objs, impulse, resting, suggest);
#ifndef NDEBUG
            std::cout << "depth: "<< ret << std::endl;
            std::cout << "t: "<< t << std::endl;
#endif
            if (ret < -depth_tolerance) {
                ode(objs, -cur_step);
                t -= cur_step;
                cur_step /= 2.0;
                flag = true;
            } else if (flag == false) {
                break;
            } else if (ret > depth_tolerance) {
                ode(objs, cur_step);
                t += cur_step;
                cur_step /= 2.0;
            } else {
                break;
            }
        }

        if (t1 - lframe > draw_step) {
            frame_count++;
            std::cout << "Frame: " << frame_count << std::endl;
            lframe = t1;
            // drawing procedure
            mon.clear();

            for (int i = 0; i < objs.size(); i++) {
                auto c = objs[i];
                Color color;
                color_map.getColorAtValue(static_cast<float>(i) / static_cast<float>(objs.size()), color);
                mon.drawPoly(objs[i].world, 1, color);
            }

            mon.refresh();
        }
    }
}
