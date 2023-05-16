#include "triangulate.h"

#include "triangle.h"

void Triangulate::triangulateSolid(RigidPoly<Scalar, Vector>& poly) {
    triangulateio input;

    input.pointlist = (double *)malloc(sizeof(double) * 2 * poly.verts.size());
    for (int i = 0; i < poly.verts.size(); i++) {
        input.pointlist[i * 2] = poly.verts[i](0);
        input.pointlist[i * 2 + 1] = poly.verts[i](1);
    }
    input.numberofpoints = poly.verts.size();
    input.pointmarkerlist = NULL;
    input.numberofpointattributes = 0;

    input.numberoftriangles = 0;
    input.numberoftriangleattributes = 0
    ;
    input.segmentlist = (int *)malloc(sizeof(int) * 2 * poly.verts.size());
    for (int i = 0; i < poly.verts.size(); i++) {
        input.segmentlist[i * 2] = i;
        input.segmentlist[i * 2 + 1] = (i + 1) % poly.verts.size();
    }
    input.numberofsegments = poly.verts.size();
    input.segmentmarkerlist = NULL;

    input.numberofholes = 0;
    input.numberofregions = 0;
    input.numberofholes = 0;

    triangulateio output;

    output.pointlist = NULL;
    output.pointattributelist = NULL;
    output.pointmarkerlist = NULL;
    
    output.trianglelist = NULL;
    output.triangleattributelist = NULL;
    output.neighborlist = NULL;

    output.segmentlist = NULL;
    output.segmentmarkerlist = NULL;

    output.edgelist = NULL;
    output.edgemarkerlist = NULL;
    output.normlist = NULL;

    char option[] = "pzqa0.05";

    triangulate(option, &input, &output, NULL);

    for (int i = 0; i < output.numberofpoints; i++) {
        TRI_REAL x, y;
        x = output.pointlist[i * 2 + 0];
        y = output.pointlist[i * 2 + 1];

        poly.inner_verts.push_back({x, y});
        poly.inner_world.push_back({x, y});
        poly.inner_world[i] += poly.center;
    }

    // triangles
    for (int i = 0; i < output.numberoftriangles; i++) {
        int v0, v1, v2;
        v0 = output.trianglelist[i * 3 + 0];
        v1 = output.trianglelist[i * 3 + 1];
        v2 = output.trianglelist[i * 3 + 2];

        poly.triangles.push_back({v0, v1, v2});
    }

    poly.bound_marks = std::vector<int>(output.pointmarkerlist, output.pointmarkerlist + output.numberofpoints);

    poly.triangulated = true;
}
