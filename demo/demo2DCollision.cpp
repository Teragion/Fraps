#include "util/color.h"
#include "util/heatmap.h"
#include "util/monitor.h"

#include "circle.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <limits>
#include <vector>

#include <Eigen/Dense>

#include <chrono>
#include <iostream>
#include <thread>

using Scalar = double;
using Vector = Vec2d;

Scalar eps = 0.7; // coefficient of restitution
Vector g = {0.0, -50.0}; // gravity
Scalar dt = 1E-3; // timestep
Scalar draw_step = 0.017; // frame update interval

Scalar tol = 1E-9; // tolerance for contact detection

Scalar tol_rel_v = 1E-2; // tolerance for contact detection
// relv > tol == moving away
// relv > -tol == resting contact
// relv < -tol == collision 

Scalar depth_tolerance = 1E-2; // tolerance for reducing timestep

using contact_type = Contact<Scalar, Vector>;

void ode(std::vector<Circle<Scalar, Vector> >& objs, Scalar dt) {
    // assuming F is computed and filled
    for (auto& o : objs) {
        if (o.fixed) {
            continue;
        }
        if (dt > 0) {
            o.m += o.F * dt;
            o.v = o.m / o.w;
            o.center += o.v * dt;
        } else {
            o.v = o.m / o.w;
            o.center += o.v * dt;
            o.m += o.F * dt;
        }
    }
}

Scalar detect_contact(std::vector<Circle<Scalar, Vector> >& objs, 
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
                if (v_rel_s > tol_rel_v) { // moving away
                    continue;
                } else if (v_rel_s > -tol_rel_v) { // resting contact
                    // compute support force
                    // if (objs[i].fixed)  {
                    //     continue;
                    // } else {
                    //     Vector N = contact.direction * dot(objs[i].F, contact.direction);
                    //     objs[i].F += N;
                    // }
                    contact.i = i;
                    contact.j = j;
                    resting.push_back(contact);
                } else { // collision
                    contact.i = i;
                    contact.j = j;
                    impulse.push_back(contact);
                }
            }
        }
    }
    return 0;
}

void process_impulse(std::vector<Circle<Scalar, Vector> >& objs, 
                     std::vector<contact_type>& impulse) {
    for (auto& contact : impulse) {
        // compute impulse exchange
        Scalar numerator = -(1 + eps) * contact.v_rel_s;
        Scalar denominator = 1 / objs[contact.i].w + 1 / objs[contact.j].w; // ignoring other terms for circle
        Scalar j = numerator / denominator;
        // J = j * direction
        objs[contact.i].m -= j * contact.direction;
        if (!objs[contact.j].fixed) {
            objs[contact.j].m += j * contact.direction;
        }
    }
}

void clear_variables(std::vector<Circle<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        clear_vec(o.v);

        clear_vec(o.F);
    }
}

void update_derived(std::vector<Circle<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        o.v = o.m / o.w;
    }
}

void apply_gravity(std::vector<Circle<Scalar, Vector> >& objs) {
    for (auto& o : objs) {
        if (o.fixed) { continue; }
        o.F += g;
    }
}

Eigen::VectorX<Scalar> fdirection(const Eigen::MatrixX<Scalar>& A, const Eigen::VectorX<int>& C, int d, int n) {
    Eigen::VectorX<Scalar> delta_f = Eigen::VectorX<Scalar>::Zero(n);
    delta_f(d) = 1;

    // construct matrix Acc
    int c = C.sum();
    Eigen::MatrixX<Scalar> Acc(c, c);
    for (int i = 0; i < n; i++) {
        int row = 0;
        if (C(i) == 0) {
            continue;
        }
        for (int j = 0; j < n; j++) {
            int col = 0;
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
    for (int i = 0; i < n; i++) {
        int row = 0;
        if (C(i) == 0) {
            continue;
        }
        Acd(row) = -A(row, d);
        row++;
    }

    // solve
    Eigen::VectorX<Scalar> x = Acc.lu().solve(Acd);

    // transfer x to delta_f
    for (int i = 0; i < n; i++) {
        if (C(i) == 1) {
            delta_f(i) = x(i);
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
#ifndef NDEBUG
        std::cout << "itr = " << itr++ << std::endl;
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

void solve_contact_forces(std::vector<Circle<Scalar, Vector> >& objs, std::vector<contact_type>& resting) {
    int n_contacts = resting.size();

    if (n_contacts == 0) { return; }

#ifndef NDEBUG
    std::cout << "n_contacts = " << n_contacts << std::endl;
#endif

    if (n_contacts == 3) {
        n_contacts = 3;
    }

    Eigen::MatrixX<Scalar> A = Eigen::MatrixX<Scalar>::Zero(n_contacts, n_contacts);
    Eigen::VectorX<Scalar> b = Eigen::VectorX<Scalar>::Zero(n_contacts);
    Eigen::VectorX<Scalar> f = Eigen::VectorX<Scalar>::Zero(n_contacts);

    // computing matrix A
    for (int i = 0; i < n_contacts; i++) {
        for (int j = 0; j < n_contacts; j++) {
            // calculating a_ij : d_i_pp's dependence on f_j
            int sign;
            Scalar mass; // mass of object of interest

            int body1 = resting[i].i;
            int body2 = resting[i].j;

            A(i, j) = 0;

            if (resting[j].i == body1) {
                // -fjnj affecting body 1
                sign = 1;
                mass = objs[body1].w;
                if (objs[body1].fixed) { sign = 0; }
                A(i, j) += sign * dot(resting[i].direction, resting[j].direction) / mass;
            } else if (resting[j].j == body1) {
                // fjnj affecting body 1
                sign = -1;
                mass = objs[body1].w;
                if (objs[body1].fixed) { sign = 0; }
                A(i, j) += sign * dot(resting[i].direction, resting[j].direction) / mass;
            } 
            if (resting[j].i == body2) {
                // -fjnj affecting body 2
                sign = -1;
                mass = objs[body2].w;
                if (objs[body2].fixed) { sign = 0; }
            } else if (resting[j].j == body2) {
                // fjnj affecting body 2
                sign = 1;
                mass = objs[body2].w;
                if (objs[body2].fixed) { sign = 0; }
            } else {
                // no effect
                continue;
            }
            A(i, j) += sign * dot(resting[i].direction, resting[j].direction) / mass;
        }
    }

    // computing vector b
    for (int i = 0; i < n_contacts; i++) {
        int body1 = resting[i].i;
        int body2 = resting[i].j;

        b(i) = -dot(resting[i].direction, (objs[body1].F / objs[body1].w) - (objs[body2].F / objs[body2].w));
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

        objs[body1].F -= f(i) * resting[i].direction;
        objs[body2].F += f(i) * resting[i].direction;
    }
}

int main() {
    // initialize objects
    std::vector<Circle<Scalar, Vector> > objs;

    // last balls are fixed
    objs.push_back(Circle<Scalar, Vector>({0, 8}, 1, 9, {0, 0}));
    objs.push_back(Circle<Scalar, Vector>({1, 5}, 2, 9, {0, 0}));

    objs.push_back(Circle<Scalar, Vector>({0, -3}, 2, 10000, {0, 0}));
    objs[2].fixed = true;
    objs.push_back(Circle<Scalar, Vector>({4, 0}, 3, 10000, {0, 0}));
    objs[3].fixed = true;
    objs.push_back(Circle<Scalar, Vector>({-4, 0}, 2, 10000, {0, 0}));
    objs[4].fixed = true;

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
        std::cout << "solve_contact_forces()"<< std::endl;
#endif
        solve_contact_forces(objs, resting);

        // fine tune the step size
        Scalar cur_step = dt / 2.0;
        ode(objs, dt);
        t += dt;
        bool flag = false;
        while (true) {
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
                mon.drawCircle(c.center(0), c.center(1), c.r, 50, color);
            }

            mon.refresh();
        }
    }
}
