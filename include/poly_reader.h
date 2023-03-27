#ifndef POLY_READER_H
#define POLY_READER_H

#include <fstream>
#include <iostream>

#include <common_defs.h>
#include <rigidpoly.h>

class PolyReader {
public:
    static RigidPoly<Scalar, Vector> readFile(const char* path) {
        std::ifstream in(path, std::ios::in);
        RigidPoly<Scalar, Vector> obj;
        if (!in) {
            std::cerr << "Cannot open " << path << std::endl;
            exit(1);
        }
        std::string line;
        while (std::getline(in, line)) {
            // check v for vertices
            if (line.substr(0, 2) == "v ") {
                std::istringstream v(line.substr(2));
                double x, y, z;
                v >> x;
                v >> y;
                v >> z;
                obj.verts.push_back({x, y, z});
            }
            // check for texture coordinate
            else if (line.substr(0, 2) == "vt") {
                // do not do any processing for now
            }
            // check for faces
            else if (line.substr(0, 2) == "f ") {
                int a, b, c;  // to store mesh index
                // int A, B, C;  // to store texture index
                // std::istringstream v;
                // v.str(line.substr(2));
                const char* chh = line.c_str();
                sscanf(chh, "f %i %i %i", &a, &b, &c);  // here it read the line start with f and store the corresponding values in the variables
                // sscanf(chh, "f %i/%i %i/%i %i/%i", &a, &A, &b, &B, &c, &C);  // here it read the line start with f and store the corresponding values in the variables

                // v>>a;v>>b;v>>c;
                a--;
                b--;
                c--;
                // A--;
                // B--;
                // C--;

                obj.faces.push_back({a, b, c});
            }
        }

        return obj;
    }

}

#endif // POLY_READER_H
