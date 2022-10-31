#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <collision_parameters.h>
#include <vector>

#include "vec.h"

/**
 * This function return the array of indices of the input 2D points on the convex hull
 * in the counterclockwise orientation.
 */
template <typename Vector>
std::vector<Vector> convex_hull2(const std::vector<Vector>& points) {
    std::vector<std::pair<int, double> > idx_cos;
    std::vector<Vector> ret;
    int p0_idx = 0;
    Vector p0 = points[0];
    for (int i = 0; i < points.size(); i++) {
        if ((p0[1] > points[i][1]) || (p0[1] == points[i][1] && p0[0] > points[i][0])) {
            p0_idx = i;
            p0 = points[i];
        }
    }

    // compute and sort points by cosine value
    Vector x_axis = {1, 0};
    for (int i = 0; i < points.size(); i++) {
        if (i == p0_idx) {
            continue;
        }
        Vector dir = points[i] - p0;
        idx_cos.emplace_back(i, x_axis.dot(dir) / norm(dir));
    }
    auto comp = [&](const std::pair<int, double> &a, const std::pair<int, double> &b) {
        if (std::abs(a.second - b.second) > alignment) {
            return a.second > b.second;
        } else {
            double dista = norm2(points[a.first] - p0);
            double distb = norm2(points[b.first] - p0);
            return dista > distb;
        }
    };
    std::sort(idx_cos.begin(), idx_cos.end(), comp);

    // check for collinear points
    for (int i = 1; i < idx_cos.size() - 1; i++) {
        if (std::abs(idx_cos[i].second - idx_cos[i + 1].second) < alignment) {
            // only keep the furthest
            idx_cos[i + 1].first = -1; // mark for delete
        }
    }
    idx_cos.emplace_back(p0_idx, 0.0);

    ret.clear();
    ret.push_back(p0);
    unsigned int stack_top = 0;
    for(auto itr = idx_cos.begin(); itr != idx_cos.end(); itr++) {
        if(itr->first == -1) { // skip if deleted
            continue;
        }
        unsigned int p3_idx = itr->first;
        while(true) {
            if(stack_top == 0) {
                break;
            }
            unsigned int p1_idx = ret[stack_top - 1], p2_idx = ret[stack_top];
            Vector p1p2 = points[p2_idx] - points[p1_idx];
            Vector p1p3 = points[p3_idx] - points[p1_idx];
            if(p1p2[0] * p1p3[1] - p1p2[1] * p1p3[0] <= 0) { // right turn or collinear
                ret.pop_back(); // pop top of the stack
                stack_top--;
            } else {
                break;
            }
        }
        ret.push_back(points[p3_idx]);
        stack_top++;
    }
}

#endif  // GEOMETRY_H
