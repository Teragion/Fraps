#include "util.h"
#include "util/color.h"
#include "util/heatmap.h"
#include "util/monitor.h"

#include "collision_parameters.h"
#include "poly_writer.h"
#include "rigidpoly.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <limits>
#include <vector>

#include <Eigen/Dense>

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

#include <cassert>

#include <png.h>

const int width = 1024;
const int height = 1024;

std::vector<Vec2d> points;
const char* path = "/Users/teragion/Models/out.poly";

Vec2d translate(double mx, double my) {
    return {(mx / (width / 2)) - 1, -(my / (height / 2) - 1)};
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) 
    {
        double xpos, ypos;
        //getting cursor position
        glfwGetCursorPos(window, &xpos, &ypos);
        std::cout << "Click at (" << xpos << " , " << ypos << ")" << std::endl;

        points.push_back(translate(xpos, ypos));
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_Z && action == GLFW_PRESS) {
        if (!points.empty()) {
            points.pop_back();
        }
    }

    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        PolyWriter::writeSolid(path, points);
    }
}

int main() {
    // initialize objects

    Monitor mon(width, height);

    Heatmap color_map;
    
    glMatrixMode(GL_PROJECTION);
    glScalef(1.0, 1.0, 1.0);

    const double t0 = glfwGetTime();
    double t = 0.0;

    ::glfwSetMouseButtonCallback(mon.window, mouse_button_callback);
    ::glfwSetKeyCallback(mon.window, key_callback);

    while (!mon.shouldClose) {
        mon.clear();

        double mxpos, mypos;

        glfwGetCursorPos(mon.window, &mxpos, &mypos);
        Vec2d sp = translate(mxpos, mypos);

        for (auto& p : points) {
            // draw circle
            mon.drawCircle(p(0), p(1), 0.002, 20);
        }

        mon.drawCircle(sp(0), sp(1), 0.002, 20);

        if (!points.empty()) {
            for (int i = 0; i < points.size() - 1; i++) {
                int j = (i + 1) % points.size();

                mon.drawLine(points[i], points[j]);
            }

            mon.drawLine(points[points.size() - 1], sp);
        }

        mon.refresh();
    }
}
