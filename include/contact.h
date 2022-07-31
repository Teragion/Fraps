#ifndef CONTACT_H
#define CONTACT_H

template <typename Scalar, typename Vector>
struct Contact {
    Vector direction;
    Scalar depth;
    Vector v_rel;

    Scalar v_rel_s;

    // indices of objects, temporary
    int i;
    int j;
};

#endif
