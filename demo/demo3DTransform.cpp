#include <arm_neon.h>
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "Eigen/src/Geometry/Quaternion.h"

#include "logger.h"
#include "obj_reader.h"
#include "rigidpoly3D.h"
#include "util/monitor.h"

#include "OpenGL/glu.h"

// #include "collision_parameters.h"

const double dt = 1E-5;

using Scalar = double;

int main() {
    Monitor mon(1024, 1024, 3);

    std::vector<RigidPoly3D<Scalar> > objs;
    objs.push_back(ObjReader::readFile("/Users/teragion/Models/box.obj"));
    // objs.push_back(ObjReader::readFile("/Users/teragion/Models/cylinder.obj"));
    // objs.push_back(ObjReader::readFile("/Users/teragion/Models/rect.obj"));
    std::cout << objs[0].I << std::endl;
    std::cout << objs[0].w << std::endl;

    for (auto& f : objs[0].verts) {
        printf("%+2.1lf %+2.1lf %+2.1lf\n", f(0), f(1), f(2));
    }

    Scalar t = 0;

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat lightpos[] = {3, 1, 2, 0.};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

    while (!mon.shouldClose) {
        mon.clear();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        // gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -100.0, 0.0, 1.0, 0.0); // default
        glTranslatef(0.0f, 0.0f, -6.0f);
        glScalef(0.5f, 0.5f, 0.5f);

        Scalar roll = 1.5707, pitch = 0, yaw = 0.707;
        Scalar scale = t;

        Eigen::Quaternion<Scalar> q;
        q = Eigen::AngleAxisd(roll * scale, Eigen::Vector3d::UnitX())
            * Eigen::AngleAxisd(pitch * scale, Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd(yaw * scale, Eigen::Vector3d::UnitZ());
        
        GLfloat cyan[] = {0.f, .8f, .8f, 1.f};
        glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
        
        for (auto& obj : objs) {
            obj.qrot = q * obj.qrot;
            obj.local_to_world();
            mon.drawPoly3D(obj);
        }

        t += dt;

        mon.refresh();
    }
}
