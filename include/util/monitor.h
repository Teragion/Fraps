#ifndef MONITOR_H
#define MONITOR_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <string>
#include <vector>

#ifdef FRAPS_3D
#include "rigidpoly3D.h"
#endif

#include "util/color.h"

class Monitor {
   public:
    Monitor();
    // Monitor(int width, int height);
    Monitor(int width, int height, int D = 2);
    ~Monitor();

    void init();
    void clear();
    void refresh();

    // 2D methods
    template<typename Vec>
    void drawLine(const Vec& v0, const Vec& v1, const float& width = 1, const Color& color = white);
    template<typename Vec>
    void drawPoly(const std::vector<Vec>& points, const float& width = 1, const Color& color = white);
    inline void drawCircle(float cx, float cy, float r, int num_segments, const Color& color = white);
    template<typename Vec>
    void drawQuad(const Vec& v0, const Vec& v1, const Vec& v2, const Vec& v3, const Color& color = white);
    template<typename Vec>
    void drawRect(const Vec& v0, const Vec& v1, const Color& color = white);
    template<typename Vec>
    void drawTri(const Vec& v0, const Vec& v1, const Vec& v2, const Color& color = white);

    // 3D methods
    #ifdef FRAPS_3D
    template<typename Scalar>
    void drawPoly3D(const RigidPoly3D<Scalar>& poly);
    template<typename Vec>
    void drawTri3D(const Vec& v0, const Vec& v1, const Vec& v2, const Color& color = white);
    #endif

    bool initialized = false;
    bool shouldClose = false;
    GLFWwindow* window = nullptr;
    double scale = 1.0;
    double trans[3] = {0, 0, 0};
    Color bgcolor = {0, 0, 0, 0};
    unsigned int width = 1024;
    unsigned int height = 1024;
    int D;
    std::string window_title = "FRAPS Monitor";
};

template<typename Vec>
void Monitor::drawLine(const Vec& v0, const Vec& v1, const float& width, const Color& color) {
    ::glLineWidth(width);
    ::glBegin(GL_LINES);
    ::glColor3f(color.r, color.g, color.b);
    ::glVertex2d(v0(0), v0(1));
    ::glVertex2d(v1(0), v1(1));
    ::glEnd();
}

template<typename Vec>
void Monitor::drawPoly(const std::vector<Vec>& points, const float& width, const Color& color) {
    ::glLineWidth(width);
    ::glBegin(GL_POLYGON);
    ::glColor3f(color.r, color.g, color.b);
    for (auto const& point : points) {
        ::glVertex2d(point(0), point(1));
    }
    ::glEnd();
}

inline void Monitor::drawCircle(float cx, float cy, float r, int num_segments, const Color& color) {
    ::glBegin(GL_POLYGON);
    ::glColor3f(color.r, color.g, color.b);
    for (int ii = 0; ii < num_segments; ii++) {
        float theta = 2.0f * M_PI * static_cast<float>(ii) / static_cast<float>(num_segments);  // get the current angle

        float x = r * cosf(theta);  // calculate the x component
        float y = r * sinf(theta);  // calculate the y component

        glVertex2f(x + cx, y + cy);  // output vertex
    }
    glEnd();
}

template<typename Vec>
void Monitor::drawQuad(const Vec& v0, const Vec& v1, const Vec& v2, const Vec& v3, const Color& color) {
    ::glBegin(GL_QUADS);
    ::glColor3f(color.r, color.g, color.b);
    ::glVertex2d(v0(0), v0(1));
    ::glVertex2d(v1(0), v1(1));
    ::glVertex2d(v2(0), v2(1));
    ::glVertex2d(v3(0), v3(1));
    ::glEnd();
}

template<typename Vec>
void Monitor::drawRect(const Vec& v0, const Vec& v1, const Color& color) {
    ::glBegin(GL_QUADS);
    ::glColor3f(color.r, color.g, color.b);
    // counterclockwise
    ::glVertex2d(v0(0), v0(1));
    ::glVertex2d(v0(0), v1(1));
    ::glVertex2d(v1(0), v1(1));
    ::glVertex2d(v1(0), v0(1));
    ::glEnd();
}

template<typename Vec>
void Monitor::drawTri(const Vec& v0, const Vec& v1, const Vec& v2, const Color& color) {
    ::glBegin(GL_TRIANGLES);
    ::glColor3f(color.r, color.g, color.b);
    ::glVertex2d(v0(0), v0(1));
    ::glVertex2d(v1(0), v1(1));
    ::glVertex2d(v2(0), v2(1));
    ::glEnd();

    // drawLine(v0, v1, 1, black);
    // drawLine(v1, v2, 1, black);
    // drawLine(v0, v2, 1, black);
}

#ifdef FRAPS_3D
template<typename Scalar>
void Monitor::drawPoly3D(const RigidPoly3D<Scalar>& poly) {
    for (int i = 0; i < poly.faces.size(); i++) {
        drawTri3D(poly.world[poly.faces[i][0]], poly.world[poly.faces[i][1]], poly.world[poly.faces[i][2]]);
    }
}

template<typename Vec>
void Monitor::drawTri3D(const Vec& v0, const Vec& v1, const Vec& v2, const Color& color) {
    ::glBegin(GL_TRIANGLES);
    ::glColor3f(color.r, color.g, color.b);
    Vec normal = (v1 - v0).cross(v2 - v0);
    ::glNormal3d(normal(0), normal(1), normal(2));
    ::glVertex3d(v0(0), v0(1), v0(2));
    ::glVertex3d(v1(0), v1(1), v1(2));
    ::glVertex3d(v2(0), v2(1), v2(2));
    ::glEnd();
}
#endif

#endif // MONITOR_H
