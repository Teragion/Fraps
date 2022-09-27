#ifndef FORCE_H
#define FORCE_H

template<typename Scalar, typename Vector>
struct Force {
    Vector dir;
    Scalar mag;
};

template<typename Scalar, typename Vector>
struct Impact {
    Vector pos;
    Vector J;
};

#endif // FORCE_H
