#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "contact.h"
#include "vec.h"

enum RIGID_TYPE { CIRCLE, POLY, TRI };

template<typename Scalar, typename Vector>
struct RigidBody {
    typedef Contact<Scalar, Vector> contact_type;

    int id;
    const RIGID_TYPE type;

    Vector center;

    // primary variables
    Scalar w; // weight
    Vector m; // momentum
    Scalar L; // angular momentum
    Scalar I; // "inertia tensor"

    // derived variables
    Vector v; // velocity (d)
    Scalar omega; // angular velocity (d)

    // computed variables
    Vector F;
    Scalar torque;

    // other parameters
    bool fixed;
};

#endif
