#ifndef RIGIDPOLY3D_H
#define RIGIDPOLY3D_H

#include <assert.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <initializer_list>
#include <limits>
#include <list>
#include <vector>

#include "Eigen/src/Geometry/Quaternion.h"
#include "contact.h"
#include "mat.h"
#include "vec.h"
// #include "section.h"

#ifndef FRAPS_USE_EIGEN
#define FRAPS_USE_EIGEN
#endif

template<typename Scalar>
struct RigidPoly3D {
#ifndef FRAPS_USE_EIGEN
    using Vector = Vec<3, Scalar>;
    using Matrix = Mat<3, 3, Scalar>;
#else
    using Vector = Eigen::Matrix<Scalar, 3, 1>;
    using Matrix = Eigen::Matrix<Scalar, 3, 3>;
#endif

    typedef Vec3i tri;
    std::vector<Vector> verts;
    std::vector<tri> faces;

#ifdef RefCode2D
    bool decomposed = false;
    typedef std::vector<Vector> Shape;
    std::vector<Shape> convex_polys;
    std::vector<Shape> convex_polys_world;
    std::vector<Vector> steinerPoints, reflexVertices;  // auxillary arrays
#endif

    // rigid body related variables
    typedef Contact<Scalar, Vector> contact_type;

    Vector center;
    Scalar rot;
    Eigen::Quaternion<Scalar> qrot = Eigen::Quaternion<Scalar>::Identity();
    Matrix rot_matrix;

    Vector center_last;
    Scalar rot_last;

    // primary variables
    Scalar w;  // weight
    Vector m;  // momentum
    Matrix L;  // angular momentum
    Matrix I;  // "inertia tensor"

    Scalar w_inv;
    Matrix I_inv;

    Scalar density = 1.0;

    // derived variables
    Vector v;      // velocity (d)
    Matrix omega;  // angular velocity (d)

    Vector v_last;
    Matrix omega_last;

    // computed variables
    Vector F;                   // force
    Matrix torque;              // torque (perpendicular to the plane)
    std::vector<Vector> world;  // verts in world coordinate

    // fracture related varibles
    // std::vector<Section> sections;
    Scalar fracture_dt = 0.01;
    Scalar stress_toleration = 2.0;
    bool to_fracture = false;
    int fracture_pos = -1;
    // Impact<Scalar, Vector> causing_impact;
    bool first_impact = true;
    bool detached = false;
    Vector new_m0;
    Vector new_m1;
    Scalar new_L0;
    Scalar new_L1;
    // std::list<Impact<Scalar, Vector> > impacts;

    // other parameters
    bool fixed = false;
    bool enable_fracture = true;
    bool deleted = false;

    RigidPoly3D() : verts(), fixed(false) {}

    RigidPoly3D(std::initializer_list<Vector> vertices, Scalar density)
        : verts(vertices), center({0.0, 0.0}), rot(0.0), m(0.0), L(0.0), density(density), fixed(false), enable_fracture(false) {
        init();
    }

    contact_type detect_collision(RigidPoly3D& other);

    void init();

    void triangulate();
    void decompose(std::vector<Vector>& poly);
    void convex_decompose();

    void correct_center();
    void init_inertia();

    // void calculate_section_mass_mmoi(Section& section);
    void init_sections();

    void apply_impulse(Vector r, Vector j, int i_s, Scalar i_r);
    bool check_fracture();
    void distribute_impulse();
    bool find_impulse_side(Scalar i_r, int i_s, int id);

    void update_rot_matrix();
    Vector local_to_world(const Vector& r);
    void local_to_world();
    Vector calculate_world_vel(Vector r);
    Vector calculate_world_vel_last(Vector r);

    void set_fixed();
};

template <typename Scalar>
void RigidPoly3D<Scalar>::init() {
    // TODO: implementation
    for (int i = 0; i < verts.size(); i++) {
        world.push_back(verts[i]);
    }

    init_inertia();
}

template<typename Vector>
inline Vector tet_mean(Vector w0, Vector w1, Vector w2, Vector w3) {
    return (w0 + w1 + w2 + w3) / 4;
}

template<typename Vector, typename Matrix>
inline Matrix translate(Matrix C, Vector x, Vector center, double mass) {
    return C + mass * (x * center.transpose() + center * x.transpose() + x * x.transpose());
}

template<typename Scalar>
void RigidPoly3D<Scalar>::init_inertia() {
    // from Johnathan Blow's tutorial 
    // "How to find the inertia tensor (or other mass properties) of a 3D solid body represented by a triangle mesh"

    Matrix C_total = Matrix::Zero();
    Scalar sum_mass = 0.0;
    Vector weighted_sum = Vector::Zero();
    Vector O = Vector::Zero();
    Matrix C_canonical;
    C_canonical << 1.0 / 60.0 , 1.0 / 120.0, 1.0 / 120.0,
                   1.0 / 120.0, 1.0 / 60.0 , 1.0 / 120.0,
                   1.0 / 120.0, 1.0 / 120.0, 1.0 / 60.0;

    Matrix A, C_prime;
    Vector w0 = Vector::Zero();
    for (int i = 0; i < faces.size(); i++) {
        // Mat A
        Vector& w1 = verts[faces[i][0]];
        Vector& w2 = verts[faces[i][1]];
        Vector& w3 = verts[faces[i][2]];
        A.col(0) = w1;
        A.col(1) = w2;
        A.col(2) = w3;

        Scalar m_tet = density * A.determinant() / 6.0;
        // C'
        C_prime = density * A.determinant() * A * C_canonical * A.transpose();

        C_total += C_prime;

        // mass
        sum_mass += m_tet;
        weighted_sum += m_tet * tet_mean(w0, w1, w2, w3);
    }

    Vector center_offset = weighted_sum / sum_mass;
    for (auto& v : verts) {
        v -= center_offset;
    }
    center = center_offset;

    w = sum_mass * density;
    w_inv = 1.0 / w;

    // translate function
    center_offset = -center_offset;
    Matrix C = translate(C_total, center_offset, center, w);
    I = C.trace() * Matrix::Identity() - C;
    I_inv = I.inverse();
}

template<typename Scalar>
void RigidPoly3D<Scalar>::update_rot_matrix() {
    rot_matrix = qrot.toRotationMatrix();
}

template<typename Scalar>
typename RigidPoly3D<Scalar>::Vector RigidPoly3D<Scalar>::local_to_world(const Vector& r) {
    return rot_matrix * r;
}

template<typename Scalar>
void RigidPoly3D<Scalar>::local_to_world() {
    update_rot_matrix();
    for (int i = 0; i < verts.size(); i++) {
        world[i] = rot_matrix * verts[i];
        // std::cout << verts[i] << std::endl;
    }
}

#endif  // RIGIDPOLY3D_H