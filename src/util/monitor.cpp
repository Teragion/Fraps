#include "util/monitor.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cstdio>
#include <iostream>
#include <string>

void gl_error_callback(int error, const char* description) {
    fprintf(stderr, "Error: %s\n", description);
}

void gl_close_callback(GLFWwindow* window) {
    Monitor* m = static_cast<Monitor*>(glfwGetWindowUserPointer(window));
    m->shouldClose = true;
}

int init_monitor(Monitor* m) {
#ifdef FLOP_DEBUG_INFO
    printf("Compiled against GLFW %i.%i.%i\n",
       GLFW_VERSION_MAJOR,
       GLFW_VERSION_MINOR,
       GLFW_VERSION_REVISION);
#endif

    if (!glfwInit()) {
        // Initialization failed
        fprintf(stderr, "Error: Cannot initialize glfw\n");
        return 1;
    }
    ::glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    ::glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    ::glfwWindowHint(GLFW_REFRESH_RATE, 60);
    ::glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);
    ::glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    // ::glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // ::glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    m->window = glfwCreateWindow(m->width, m->height, m->window_title.c_str(), NULL, NULL);
    if (m->window == NULL) {
        fprintf(stderr, "Error: Cannot create window\n");
        return 1;
    }

    // callbacks
    ::glfwSetErrorCallback(gl_error_callback);
    ::glfwSetWindowCloseCallback(m->window, gl_close_callback);

    ::glfwMakeContextCurrent(m->window);
    ::gladLoadGL();
    ::glfwSwapInterval(1);

    m->initialized = true;

    ::glMatrixMode(GL_PROJECTION);
    // ::glLoadIdentity();
    // ::glOrtho(0, m->width, m->height, 0, 0, 1);

#ifdef FLOP_DEBUG_INFO
    std::cout << "OpenGL version: " << reinterpret_cast<const char *>(glGetString(GL_VERSION)) << std::endl;
#endif

    return 0;
}

void dest_monitor(Monitor* m) {
    ::glfwTerminate();
}

Monitor::Monitor() {
    Monitor(1024, 1024);
}

Monitor::Monitor(int width, int height) : width(width), height(height) {
    init();
}

Monitor::~Monitor() {
    fprintf(stdout, "Closing.\n");
    dest_monitor(this);
}

void Monitor::init() {
    if (init_monitor(this)) {
        printf("Error.\n");
        return;
    }

    // allow reflecting from callbacks
    ::glfwSetWindowUserPointer(window, this);
}

void Monitor::clear() {
    ::glClearColor(
      static_cast<GLclampf>(bgcolor.r),
      static_cast<GLclampf>(bgcolor.g),
      static_cast<GLclampf>(bgcolor.b),
      static_cast<GLclampf>(bgcolor.a));
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Monitor::refresh() {
    ::glFlush();
    ::glfwSwapBuffers(window);
    ::glfwPollEvents();
}
