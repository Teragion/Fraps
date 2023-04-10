#ifndef POLY_READER_H
#define POLY_READER_H

#include <fstream>
#include <iostream>

#include <common_defs.h>
#include <rigidpoly.h>

class PolyReader {
public:
    static RigidPoly<Scalar, Vector> readSolid(const char* path) {
        std::ifstream in(path, std::ios::in);
        RigidPoly<Scalar, Vector> obj;
        if (!in) {
            std::cerr << "Cannot open " << path << std::endl;
            exit(1);
        }
        std::string line;
        std::getline(in, line);

        std::istringstream v(line);
        int verts;
        v >> verts;

        for (int i = 0; i < verts; i++) {
            std::getline(in, line);
            std::istringstream v(line);

            int id;
            double x, y;

            v >> id;
            v >> x;
            v >> y;

            obj.verts.push_back({x, y});
        }

        in.close();

        obj.init();

        return obj;
    }

};

#endif // POLY_READER_H
