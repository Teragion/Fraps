#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H

#include "vec.h"
#define FRAPS_2D
// #define FRAPS_3D

#ifdef FRAPS_3D
    #ifndef FRAPS_USE_EIGEN
    #define FRAPS_USE_EIGEN
    #endif
#endif

using Scalar = double;
#ifdef FRAPS_USE_EIGEN
#include <Eigen/Dense>
using Vector = Eigen::Matrix<double, 3, 1>;
#else
using Vector = Vec2d;
#endif

#endif // COMMON_DEFS_H
