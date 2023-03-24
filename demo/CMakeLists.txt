
list(TRANSFORM GLAD_SRCS PREPEND ../)

find_package(PNG REQUIRED)

add_executable(Fracture2DDemo
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    demo2DFracture.cpp
)

target_link_libraries(Fracture2DDemo
  OpenGL::GL
  glfw
  # ${MKL_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${PNG_LIBRARY}
  # OpenMP::OpenMP_CXX
  # Eigen3::Eigen
)

add_executable(Collision2DDemo
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    demo2DCollision.cpp
)

target_link_libraries(Collision2DDemo
  OpenGL::GL
  glfw
  # ${MKL_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${PNG_LIBRARY}
  # OpenMP::OpenMP_CXX
  # Eigen3::Eigen
)

add_executable(PolyCollision2DDemo
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    demo2DPolyCollision.cpp
)

target_link_libraries(PolyCollision2DDemo
  OpenGL::GL
  glfw
  # ${MKL_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${PNG_LIBRARY}
  # OpenMP::OpenMP_CXX
  # Eigen3::Eigen
)

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
)

add_executable(Transform3DDemo
    ${GLAD_SRCS}
    ${FRAPS3D_SRCS}
    ${FRAPS3D_HDRS}
    demo3DTransform.cpp
)

target_compile_options(Transform3DDemo 
                       PUBLIC -DFRAPS_USE_EIGEN)

set_target_properties(Transform3DDemo PROPERTIES
                      VERIFY_INTERFACE_HEADER_SETS TRUE
)

target_link_libraries(Transform3DDemo
    OpenGL::GL
    glfw
    # ${MKL_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${PNG_LIBRARY}
    # OpenMP::OpenMP_CXX
    # Eigen3::Eigen
)


