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

    int i_s;
    Scalar i_r;
};

#endif // FORCE_H
