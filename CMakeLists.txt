cmake_minimum_required(VERSION 3.12)
project(Fraps)

enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_CXX_EXTENSIONS OFF)
if(MSVC)
  # GTest is static compiled by default
  # set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")  
  set(CMAKE_CXX_FLAGS "/W4")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
elseif(APPLE)
  # somehow apple does not search this directory
  # include_directories(/usr/local/include) 
  # add_link_options("-L/usr/local/lib")
  add_link_options("-L/opt/homebrew/opt/llvm/lib")
  add_compile_options("-I/opt/homebrew/opt/llvm/include")
endif()
add_compile_options(-Wall)
add_compile_options(-xc++)

# Coloring
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color=always)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
endif()

################################
# release flags 
message(CMAKE_CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -msse3 -mfpmath=sse -mavx -fopenmp ")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -msse3 -mfpmath=sse -mavx ")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -fopenmp=libomp ")
  # set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -msse3 -mfpmath=sse -mavx -fopenmp=libomp ")
elseif(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "/MD /O2 /Ob2 /DNDEBUG /arch:AVX2 /fp:fast ")
endif()
message(CMAKE_CXX_FLAGS_RELEASE="${CMAKE_CXX_FLAGS_RELEASE}")

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_compile_definitions(FLOP_DEBUG_INFO=1)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_compile_definitions(FLOP_DEBUG_INFO=1)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# OpenMP
find_package(OpenMP REQUIRED)

set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# include directories 
include_directories(include/)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Eigen
include_directories(${CMAKE_SOURCE_DIR}/imported/eigen/)
find_package(Eigen3 REQUIRED)

# OpenGL
find_package(OpenGL REQUIRED)

# libMesh
find_package(LIBMESH REQUIRED)

# Intel MKL
# Unofficial
# find_package(MKL REQUIRED)
# include_directories(${MKL_INCLUDE_DIRS})
# message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
# Official
set(MKL_THREADING "intel_thread")
set(MKL_INTERFACE "lp64")
find_package(MKL CONFIG)
message(STATUS "${MKL_IMPORTED_TARGETS}")
if(${MKL_FOUND})
  set(MKL_LINK_OPTION "MKL::MKL")
endif()
include_directories(${MKL_INCLUDE})


# glad 
set(GLAD_INCLUDE_DIR "imported")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)
include_directories(${GLAD_INCLUDE_DIR})
set(GLAD_SRCS "imported/glad/glad.h" "imported/src/glad.c")

# glfw
find_package(glfw3 3.3 REQUIRED)

# check includes 
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

file(GLOB_RECURSE FRAPS_SRCS src/*.cpp)
file(GLOB_RECURSE FRAPS_HDRS include/*.h)

set(FRAPS3D_HDRS
  include/util/heatmap.h
  include/util/monitor.h
  
  include/array1.h
  include/array2_utils.h
  include/array2.h
  include/circle.h
  include/collision_parameters.h
  include/contact.h
  include/epa.h
  include/force.h
  include/fracture.h
  include/geometry.h
  include/logger.h
  include/mat.h
  include/math_kernel.h
  include/obj_reader.h
  include/rigidpoly3D.h
  include/util.h
  include/vec.h
)

set(FRAPS3D_SRCS
  src/util/heatmap.cpp
  src/util/monitor.cpp
)

list(TRANSFORM FRAPS3D_SRCS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
list(TRANSFORM FRAPS3D_HDRS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

# Testing
set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
include(CTest)
add_subdirectory(imported/googletest)

# Triangle
option(BUILD_EXECUTABLE "Build triangle executable" OFF)
option(BUILD_LIBRARY    "Build libtriangle library" ON)
add_subdirectory(imported/triangle)
include_directories(imported/triangle)

add_subdirectory(tests)

add_subdirectory(demo)
add_subdirectory(tool)
