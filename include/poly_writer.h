#ifndef POLY_WRITER_H
#define POLY_WRITER_H

#include <fstream>
#include <vector>

#include <vec.h>

// following https://www.cs.cmu.edu/~quake/triangle.poly.html
class PolyWriter {
public:
    static void writeSolid(const char* path, std::vector<Vec2d> points) {
        std::ofstream out;
        out.open(path, std::ios::out);
        // First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
        out << points.size() << " 2 0 0\n";
        for (int i = 0; i < points.size(); i++) {
            out << i + 1 << " " << points[i](0) << " " << points[i](1) << "\n";
        }

        for (int i = 0; i < points.size(); i++) {
            int j = (i + 1) % points.size();
            out << i + 1 << " " << j + 1 << "\n";
        }

        out << "0\n"; // holes

        out.close();
    }
};

#endif // POLY_WRITER_H