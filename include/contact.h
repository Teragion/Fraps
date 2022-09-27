#ifndef CONTACT_H
#define CONTACT_H

template <typename Scalar, typename Vector>
struct Contact {
    Vector p;
    Vector direction;
    Scalar depth;
    Vector v_rel;

    Vector r_i;
    Vector r_j;

    Scalar v_rel_s;

    bool vert_at_i;
    bool f2f;
    Vector f2f_offset;

    // indices of objects, temporary
    int i;
    int j;

    bool processed = false;
};

#endif // CONTACT_H
