#ifndef STRESS_H
#define STRESS_H

#include <iostream>

#include "libmesh/point.h"
#include "rigidpoly.h"
#include "section.h"
#include "util/heatmap.h"
#include "util/monitor.h"
#include "vec.h"

template<typename StressFunc>
double calculate_point_stress(StressFunc& sigma00, StressFunc& sigma01, StressFunc& sigma11, const Vec2d& p, const Vec2d& dir) {
    auto s00 = sigma00(libMesh::Point(p(0), p(1), 0));
    auto s01 = sigma01(libMesh::Point(p(0), p(1), 0));
    auto s11 = sigma11(libMesh::Point(p(0), p(1), 0));

    const Vec2d stress{dir(0) * s00 + dir(1) * s01, dir(0) * s01 + dir(1) * s11};
    return norm(stress);
}

template<typename StressFunc>
double max_segment_stress(StressFunc& sigma00, StressFunc& sigma01, StressFunc& sigma11, const Section& s) {
    constexpr int sampling_points = 50;

    Vec2d l = s.j0.p;
    Vec2d r = s.j1.p;

    double max_stress = -1;

    for (int i = 0; i < sampling_points; i++) {
        Vec2d p = lerp(l, r, static_cast<double>(i) / (sampling_points - 1));
        double stress = calculate_point_stress(sigma00, sigma01, sigma11, p, s.normal);
        if (stress > max_stress) {
            max_stress = stress;
        }
    }

    return max_stress;
}

template<typename StressFunc, typename Scalar, typename Vector>
void calculate_sections_stress(StressFunc& sigma00, StressFunc& sigma01, StressFunc& sigma11, RigidPoly<Scalar, Vector>& obj) {
    for (int i = 0; i < obj.sections.size(); i++) {
        double stress = max_segment_stress(sigma00, sigma01, sigma11, obj.sections[i]);
        double ratio = stress / obj.stress_toleration;
        obj.sections[i].ratio = std::max(ratio, obj.sections[i].ratio);
    }
}

template<typename Scalar, typename Vector>
void visualize_stress(RigidPoly<Scalar, Vector>& obj) {
    Monitor mon(1024, 1024);

    Heatmap color_map;

    glMatrixMode(GL_PROJECTION);
    glScalef(1.0, 1.0, 1.0);

    uint8_t* pixels = new uint8_t[1024 * 1024 * 3];

    unsigned int frame_count = 0;

    const double t0 = glfwGetTime();
    double lframe = t0;
    double t = 0.0;
    Scalar suggest;

    while (!mon.shouldClose) {
        const double t1 = glfwGetTime();

        if (t1 - lframe > draw_step) {
            frame_count++;
            //            std::cout << "Frame: " << frame_count << std::endl;
            lframe = t1;
            // drawing procedure
            mon.clear();

            if (obj.triangulated) {
                // for (int j = 0; j < obj.convex_polys_world.size(); j++) {
                //     mon.drawPoly(obj.convex_polys_world[j], 1, white);
                // }

                for (int j = 0; j < obj.triangles.size(); j++) {
                    const auto& tri = obj.triangles[j];
                    mon.drawTri(obj.inner_world[tri(0)], obj.inner_world[tri(1)], obj.inner_world[tri(2)], white);
                }
            } else {
                mon.drawPoly(obj.world, 1, white);
            }

            Scalar max_r;
            for (int i = 0; i < obj.sections.size(); i++) {
                if ((obj.sections[i].ratio) > max_r) {
                    max_r = obj.sections[i].ratio;
                }
            }

            for (int i = 0; i < obj.sections.size(); i++) {
                Color color;
                color_map.getColorAtValue(static_cast<float>(obj.sections[i].ratio) / max_r, color);
                mon.drawLine(obj.center + obj.sections[i].j0.p, obj.center + obj.sections[i].j1.p, 1, color);
            }

            mon.refresh();
        }
    }
}

#endif
