#ifndef POLYGON_H
#define POLYGON_H

#include <initializer_list>
#include <vector>

template<typename Scalar, typename Vector>
struct Polygon {
    std::vector<Vector> verts;

    Polygon() 
        : verts() {}
    
    Polygon(std::initializer_list<Vector> init_list)
        : verts(init_list) {}
};

#endif
