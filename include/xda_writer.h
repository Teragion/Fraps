#ifndef POLY_WRITER_H
#define POLY_WRITER_H

#include <fstream>
#include <vector>

#include <vec.h>
#include <rigidpoly.h>

class XDAWriter {
public:
    static void writeTriangleMesh(const char* path, const RigidPoly<Scalar, Vector>& poly) {
        std::ofstream out;
        out.open(path, std::ios::out);

        // metadata
        out << "libMesh-0.7.0+\n";
        out << poly.triangles.size() << " # Num. Elements\n";
        out << poly.inner_verts.size() << " # Num. Nodes\n";

        out << "." << "# boundary condition specification file\n";
        out << "n/a" << "# subdomain id specification file\n";
        out << "n/a" << "# processor id specification file\n";
        out << "n/a" << "# p-level specification file\n";

        out << poly.triangles.size() << " # n_elem at level 0, [ type (n0 ... nN-1) ]\n";
        for (int i = 0; i <  poly.triangles.size(); i++) {
            out << 3 << " " << poly.triangles[i](0) << " " << poly.triangles[i](1) << " " << poly.triangles[i](2) << "\n";
        }

        for (int i = 0; i < poly.inner_verts.size(); i++) {
            out << poly.inner_verts[i](0) << " " << poly.inner_verts[i](1)<< "\n";
        }

        out << 0 << " # Num. Boundary Conds.\n";

        out.close();
    }
};

#endif // XDA_WRITER_H