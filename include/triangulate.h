#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include "rigidpoly.h"

// A wrapper based on libigl/Triangle

class Triangulate {
public:
    static void triangulateSolid(RigidPoly<Scalar, Vector>& poly);
};


#endif // TRIANGULATE_H
