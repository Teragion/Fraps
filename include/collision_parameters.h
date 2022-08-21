#ifndef COLLISION_PARAMETERS_H
#define COLLISION_PARAMETERS_H

#include "vec.h"

using Scalar = double;
using Vector = Vec2d;

const Scalar eps = 0.2; // coefficient of restitution
const Vector g = {0.0, -10.0}; // gravity
const Scalar dt = 5E-5; // timestep
const Scalar draw_step = 0.017; // frame update interval

const Scalar tol = 1E-6; // tolerance for contact detection

const Scalar tol_rel_v = 1E-4; // tolerance for contact detection
// relv > tol == moving away
// relv > -tol == resting contact
// relv < -tol == collision 

const Scalar depth_tolerance = 1E-1; // tolerance for reducing timestep

#endif
