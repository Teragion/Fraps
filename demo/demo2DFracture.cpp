#include "util/color.h"
#include "util/heatmap.h"
#include "util/monitor.h"

#include "force.h"
#include "polygon.h"
#include "vec.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <random>

#include <png.h>

constexpr double eps = 0.000001;

struct Intersection {
    Vec<2, double> p;
    unsigned int idx_v0;
    unsigned int idx_v1;
};

struct Section {
    Intersection j0;
    Intersection j1;
    Vec<2, double> p;
    Vec<2, double> normal;
    double length;
};

/**
 * @brief calculate moment of inertia of a cross section length l
 * 
 * @param l
 * @return double 
 */
double moment_of_inertia(double l) {
    double r = l / 2;
    return 2 * (r * r * r / 3);
}

/**
 * @brief finds intersection between a line (described by normal and offset) and
 * a line segment
 * 
 * @param normal 
 * @param offset 
 * @param v1 
 * @param v2 
 * @return Vec<2, double> 
 */
bool intersection_normal_segment(const Vec<2, double>& normal, const double& offset, 
                                 const Vec<2, double>& v0, const Vec<2, double>& v1, Vec<2, double>& ret) {
    bool p0 = dot(v0, normal) < offset;
    bool p1 = dot(v1, normal) < offset;

    if (p0 == p1) { return false; }

    Vec<2, double> dir = v0 - v1;
    Vec<2, double> normal2 = {dir[1], -dir[0]};
    double offset2 = dot(normal2, v0);

    ret = {(offset * normal2[1] - normal[1] * offset2) / (normal[0] * normal2[1] - normal[1] * normal2[0]), 
           (offset * normal2[0] - normal[0] * offset2) / (normal[1] * normal2[0] - normal[0] * normal2[1])};

    return true;
}

inline Vec<2, double> find_edge_normal(unsigned int idx_v0, unsigned int idx_v1, const Polygon<double, Vec<2, double> >& shape) {
    Vec<2, double> ret{-(shape.verts[idx_v0][1] - shape.verts[idx_v1][1]), shape.verts[idx_v0][0] - shape.verts[idx_v1][0]};
    normalize(ret);
    return ret;
}

Vec<2, double> find_advance_pd(const Vec<2, double>& n, const Vec<2, double>& n0) {
    double det_inv = 1.0 / (n0[0] * n[1] - n0[1] * n[0]);
    return {det_inv * (-n0[1]), det_inv * n0[0]};
}

bool construct_cross_section(const Intersection& i0, const Intersection& i1, const Polygon<double, Vec<2, double> >& shape, Section& ret) {
    // find dpde
    Vec<2, double> n{-(i0.p[1] - i1.p[1]), i0.p[0] - i1.p[0]}; // normal of line of raw intersection
    Vec<2, double> n0 = find_edge_normal(i0.idx_v0, i0.idx_v1, shape);
    Vec<2, double> n1 = find_edge_normal(i1.idx_v0, i1.idx_v1, shape);

    if (dot(n, n0) < 0) { n0 = -n0; }
    if (dot(n, n1) < 0) { n1 = -n1; }

    if (std::abs(n0[0] * n[1] - n0[1] * n[0]) < eps) { return false; }
    Vec<2, double> dqde = find_advance_pd(n, n0);

    if (std::abs(n1[0] * n[1] - n1[1] * n[0]) < eps) { return false; }
    Vec<2, double> drde = find_advance_pd(n, n1);

    Vec<2, double> dpde = (dqde + drde) * 0.5;

    Vec<2, double> p = (i0.p + i1.p) * 0.5;

    std::cout << "dpde " << dpde[0] << " " << dpde[1] << std::endl;
    std::cout << "p " << p[0] << " " << p[1] << std::endl;

    ret.normal = dpde;

    // reconstruct intersections
    ret.j0.idx_v0 = i0.idx_v0;
    ret.j0.idx_v1 = i0.idx_v1;
    if (!intersection_normal_segment(dpde, dot(dpde, p), shape.verts[i0.idx_v0], shape.verts[i0.idx_v1], ret.j0.p)) { return false; }

    ret.j1.idx_v0 = i1.idx_v0;
    ret.j1.idx_v1 = i1.idx_v1;
    if (!intersection_normal_segment(dpde, dot(dpde, p), shape.verts[i1.idx_v0], shape.verts[i1.idx_v1], ret.j1.p)) { return false; }

    ret.length = dist(ret.j0.p, ret.j1.p);
    ret.p = p;

    std::cout << "section found" << std::endl;
    std::cout << "v0" << ret.j0.p[0] << " " << ret.j0.p[1] << std::endl;
    std::cout << "v1" << ret.j1.p[0] << " " << ret.j1.p[1] << std::endl;

    return true;
}

std::vector<Section> find_cross_sections(const Polygon<double, Vec<2, double> >& shape) {
    std::vector<Section> ret;
    const int num_directions = 12;
    const int num_offsets = 10;
    const double left = -2.0;
    const double right = 2.0;

    for (int i = 0; i < num_directions; i++) {
        const double theta = i * (M_PI / static_cast<double>(num_directions));
        for (int j = 0; j < num_offsets; j++) {
            const double offset = lerp(left, right, j / static_cast<double>(num_offsets));

            const Vec<2, double> normal{std::cos(theta), std::sin(theta)};

            std::vector<Intersection> intersections;

            for (unsigned int k = 0; k < shape.verts.size(); k++) {
                Vec<2, double> p;
                if (!intersection_normal_segment(normal, offset, shape.verts[k], 
                                                 shape.verts[(k + 1) % shape.verts.size()], p)) {
                    continue;
                }

                intersections.push_back({p, k, static_cast<unsigned int>((k + 1) % shape.verts.size())});
            }
            std::cout << std::endl << normal[0] << " " << normal[1] << " " << offset << std::endl;
            for (auto const & i : intersections) {
                std::cout << i.p[0] << " " << i.p[1] << std::endl;
            }

            assert(intersections.size() % 2 == 0);

            std::sort(intersections.begin(), intersections.end(), [](const Intersection& a, const Intersection& b) {
                if (std::abs(a.p[0] - b.p[0]) < eps) { return a.p[1] > b.p[1]; }
                else { return a.p[0] > b.p[0]; }
            });

            for (unsigned int k = 0; k < intersections.size(); k += 2) {
                Section s;
                if (construct_cross_section(intersections[k], intersections[k + 1], shape, s)) {
                    ret.push_back(s);
                }
            }
        }
    }

    return ret;
}

int main() {
    // init polygon
    Polygon<double, Vec<2, double> > H_shape{
        {-1.5, -1.5},
        {-0.5, -1.5},
        {-0.5, -0.5},
        {0.5, -0.5},
        {0.5, -1.5},
        {1.5, -1.5},
        {1.5, 1.5},
        {0.5, 1.5},
        {0.5, 0.5},
        {-0.5, 0.5},
        {-0.5, 1.5},
        {-1.5, 1.5}};

    double perimeter = 0.0;
    for (int i = 0; i < H_shape.verts.size(); i++) {
        perimeter += dist(H_shape.verts[i], H_shape.verts[(i + 1) % H_shape.verts.size()]);
    }

    std::random_device rd; 
    std::mt19937 gen(rd()); 
    // gen.seed(0);
    std::uniform_real_distribution<> dis01(0, 1);

    Monitor mon(1024, 1024);

    Heatmap stress_map;

    std::vector<Section> cross_sections = find_cross_sections(H_shape);

    const double t0 = glfwGetTime();
    
    glMatrixMode(GL_PROJECTION);
    glScalef(0.5, 0.5, 0);

    uint8_t *pixels = new uint8_t[1024 * 1024 * 3];

    unsigned int frame_count = 0;

    while (!mon.shouldClose) {
        frame_count++;

        double t1 = glfwGetTime();
        // generate impact position
        double t = (t1 - t0) / 10.0;
        t = t - std::floor(t);
        int seg_idx = 0;
        t *= perimeter; 


        for (int i = 0; i < H_shape.verts.size(); i++) {
            double segment = dist(H_shape.verts[i], H_shape.verts[(i + 1) % H_shape.verts.size()]);
            if (t >= segment) {
                t -= segment;
            } else {
                t = t / segment;
                seg_idx = i;
                break;
            }
        }

        Vec<2, double> pos = lerp(H_shape.verts[seg_idx], H_shape.verts[(seg_idx + 1) % H_shape.verts.size()], t);

        // generate impact force
        double theta = (t1 - t0) / 4.0;
        theta = M_PI * (std::sin(theta));
        Vec<2, double> dir{std::cos(theta), std::sin(theta)};
        double mag = std::sin((t1 - t0) / 3.5);
        Impact<double, Vec<2, double> > impact{pos, {dir, mag}};

        mon.clear();

        mon.drawPoly(H_shape.verts);

        for (auto const &s : cross_sections) {
            Vec<2, double> r = impact.pos - s.p;
            double ext_moment = std::abs(cross(r, impact.F.dir) * impact.F.mag);
            double max_stress = ext_moment * (s.length / 2) / moment_of_inertia(s.length);
            Color seg_color;
            stress_map.getColorAtValue(max_stress / 10, seg_color);

            mon.drawLine(s.j0.p, s.j1.p, 1, seg_color);
        }

        mon.drawLine(pos, pos + dir * mag);
        
        // write to png
        glReadPixels(0, 0, 1024, 1024, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);
        for (int j = 0; j * 2 < 1024; ++j) {
            int x = j * 1024 * 3;
            int y = (1024 - 1 - j) * 1024 * 3;
            for (int i = 1024 * 3; i > 0; --i) {
                uint8_t tmp = pixels[x];
                pixels[x] = pixels[y];
                pixels[y] = tmp;
                ++x;
                ++y;
            }
        }
        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        png_infop info = png_create_info_struct(png);
        if (!info) {
            png_destroy_write_struct(&png, &info);
            return -1;
        }
        std::string s = "output/Frame" + std::to_string(frame_count) + ".png";
        FILE *fp = fopen(s.c_str(), "wb");
        if (!fp) {
            png_destroy_write_struct(&png, &info);
            return -1;
        }

        png_init_io(png, fp);
        png_set_IHDR(png, info, 1024, 1024, 8 , PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
        if (!palette) {
            fclose(fp);
            png_destroy_write_struct(&png, &info);
            return -1;
        }
        png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
        png_write_info(png, info);
        png_set_packing(png);

        png_bytepp rows = (png_bytepp)png_malloc(png, 1024 * sizeof(png_bytep));
        for (int i = 0; i < 1024; ++i)
            rows[i] = (png_bytep)(pixels + (1024 - i - 1) * 1024 * 3);

        png_write_image(png, rows);
        png_write_end(png, info);
        png_free(png, palette);
        png_destroy_write_struct(&png, &info);

        fclose(fp);
        delete[] rows;
        
        mon.refresh();
    }

}
