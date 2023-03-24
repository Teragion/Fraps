#ifndef AABB_BOUNDING_H
#define AABB_BOUNDING_H

#include "common_defs.h"

struct AABBBounding {
    Vector vmin, vmax;

    Scalar xa() { return vmin.x(); }
    Scalar ya() { return vmin.y(); }
    Scalar za() { return vmin.z(); }

    Scalar xb() { return vmax.x(); }
    Scalar yb() { return vmax.y(); }
    Scalar zb() { return vmax.z(); }

};

#endif // AABB_BOUDNING_H
