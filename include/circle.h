#ifndef CIRCLE_H
#define CIRCLE_H

#include "contact.h"
#include "vec.h"

template<typename Scalar, typename Vector>
struct Circle {
    typedef Contact<Scalar, Vector> contact_type;

    Vector center;
    Scalar r;

    // primary variables
    Scalar w; // weight
    Vector m; // momentum

    // derived variables
    Vector v; // velocity (d)

    // computed variables
    Vector F;

    // other parameters
    bool fixed;

    Circle() 
        : center(), r(0.0) {}
    
    Circle(Vector center, Scalar r)
        : center(center), r(r) {}

    Circle(Vector center, Scalar r, Scalar w, Vector m)
        : center(center),
        r(r),
        w(w),
        m(m) {}

    // detect collision
    contact_type detect_collision(Circle& other);

    // update derived variables
    void update();
};

template<typename Scalar, typename Vector>
typename Circle<Scalar, Vector>::contact_type 
Circle<Scalar, Vector>::detect_collision(Circle<Scalar, Vector>& other) {
    Circle<Scalar, Vector>::contact_type ret;
    ret.depth = norm(other.center - center) - other.r - r;
    ret.direction = normalized(other.center - center);
    ret.v_rel = other.v - v;
    return ret;
}

template<typename Scalar, typename Vector>
void Circle<Scalar, Vector>::update() {
    // velocity
    v = m / w;
}

#endif // CIRCLE_H
