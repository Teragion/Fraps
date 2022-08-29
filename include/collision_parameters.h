#ifndef COLLISION_PARAMETERS_H
#define COLLISION_PARAMETERS_H

#include "vec.h"

using Scalar = double;
using Vector = Vec2d;

const Scalar eps = 0.4; // coefficient of restitution
const Vector g = {0.0, -10.0}; // gravity
const Scalar dt = 1E-4; // timestep
const Scalar draw_step = 0.017; // frame update interval

const Scalar tol = 1E-6; // tolerance for contact detection

const Scalar tol_rel_v = 1E-2; // tolerance for contact detection
// relv > tol == moving away
// relv > -tol == resting contact
// relv < -tol == collision 

const Scalar depth_tolerance = 5E-3; // tolerance for reducing timestep

#endif