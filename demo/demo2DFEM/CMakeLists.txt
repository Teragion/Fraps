
list(TRANSFORM GLAD_SRCS PREPEND ../../)

find_package(PNG REQUIRED)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(LIBMESH_CXXFLAGS )

add_executable(FEM2DDemo
    ${GLAD_SRCS}
    ${FRAPS_SRCS}
    ${FRAPS_HDRS}
    elasticity_system.h
    elasticity_system.cpp
    demo2DFEM.cpp
)

target_include_directories(FEM2DDemo
    PUBLIC .
    PUBLIC ${LIBMESH_INCLUDE_DIR}
)

target_compile_options(FEM2DDemo PUBLIC
    "-felide-constructors" "-Qunused-arguments" "-Wunused-parameter" "-Wunused" "-DNDEBUG")

target_link_options(FEM2DDemo PUBLIC
    "-Wl,-rpath,/Users/teragion/libmesh/lib" "-L/Users/teragion/libmesh/lib" "-lmesh_opt" "-ltimpi_opt" "-lz")

target_link_libraries(FEM2DDemo
    OpenGL::GL
    glfw
    # ${MKL_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${PNG_LIBRARY}
    # OpenMP::OpenMP_CXX
    # Eigen3::Eigen
    ${LIBMESH_LIBRARIES}
    triangle
    MPI::MPI_CXX
)
