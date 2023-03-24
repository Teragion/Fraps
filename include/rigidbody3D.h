#ifndef RIGIDBODY3D_H
#define RIGIDBODY3D_H

#include "contact.h"
#include "mat.h"
#include "vec.h"

template<typename Scalar>
struct RigidBody3D {
    using Vector = Vec<3, Scalar>;
    using Matrix = Mat<3, 3, Scalar>;

    typedef Contact<Scalar, Vector> contact_type;

    int id;

    Vector center;

    // primary variables
    Scalar w;  // weight
    Vector m;  // momentum
    Scalar L;  // angular momentum
    Matrix I;  // "inertia tensor"

    // derived variables
    Vector v;      // velocity (d)
    Scalar omega;  // angular velocity (d)

    // computed variables
    Vector F;
    Scalar torque;

    // other parameters
    bool fixed;
};

#endif  // RIGIDBODY3D_H
