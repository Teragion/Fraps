#ifndef MONITOR_H
#define MONITOR_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <string>

#include "util/color.h"

class Monitor {
public:
    Monitor();
    Monitor(int width, int height);
    ~Monitor();

    void init();
    void clear();
    void refresh();

    template<typename Vec> void drawLine(const Vec& v0, const Vec& v1, const float& width = 1, const Color& color = white);
    template<typename Vec> void drawPoly();
    template<typename Vec> void drawQuad(const Vec& v0, const Vec& v1, const Vec& v2, const Vec& v3, const Color& color = white);
    template<typename Vec> void drawRect(const Vec& v0, const Vec& v1, const Color& color = white);
    template<typename Vec> void drawTri(const Vec& v0, const Vec& v1, const Vec& v2, const Color& color = white);

    bool initialized = false;
    bool shouldClose = false;
    GLFWwindow *window = nullptr;
    double scale = 1.0;
    double trans[3] = {0, 0, 0};
    Color bgcolor = {0, 0, 0, 0};
    unsigned int width = 1024;
    unsigned int height = 1024;
    std::string window_title = "FLOP Monitor";
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
void Monitor::drawPoly() {
    // TODO: implement
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
}

#endif
