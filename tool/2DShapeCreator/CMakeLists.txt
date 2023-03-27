
list(TRANSFORM GLAD_SRCS PREPEND ../)

add_executable(2DShapeCreator
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    2DShapeCreator.cpp
)

target_link_libraries(2DShapeCreator
    OpenGL::GL
    glfw
    # ${MKL_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${PNG_LIBRARY}
    # OpenMP::OpenMP_CXX
    # Eigen3::Eigen
    triangle
)


