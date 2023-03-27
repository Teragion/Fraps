
list(TRANSFORM GLAD_SRCS PREPEND ../../)

find_package(PNG REQUIRED)

add_executable(RigidFracture2DDemo
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    demo2DRigidFracture.cpp
)

target_link_libraries(RigidFracture2DDemo
    OpenGL::GL
    glfw
    # ${MKL_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${PNG_LIBRARY}
    # OpenMP::OpenMP_CXX
    # Eigen3::Eigen
    triangle
)


