// C++ include files that we need
#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

// libMesh includes
#include "elasticity_system.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/diff_solver.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/eigen_sparse_linear_solver.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_function.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/newmark_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/perf_log.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/steady_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"
#include "poly_reader.h"
#include "rigidpoly.h"
#include "triangulate.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Matrix and right-hand side assemble
void assemble_elasticity(EquationSystems& es, const std::string& system_name);

// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
Real eval_elasticity_tensor(unsigned int i, unsigned int j, unsigned int k, unsigned int l);

void compute_stresses(EquationSystems& es);

void write_xda(const RigidPoly<Scalar, Vector>& poly) {
    std::string path = "/Users/teragion/Models/out.xda";

    std::ofstream out;
    out.open(path, std::ios::out);

    // metadata
    out << "libMesh-0.7.0+\n";
    out << poly.triangles.size() << " # Num. Elements\n";
    out << poly.inner_verts.size() << " # Num. Nodes\n";

    out << "."
        << "# boundary condition specification file\n";
    out << "n/a"
        << "# subdomain id specification file\n";
    out << "n/a"
        << "# processor id specification file\n";
    out << "n/a"
        << "# p-level specification file\n";

    out << poly.triangles.size() << " # n_elem at level 0, [ type (n0 ... nN-1) ]\n";
    for (int i = 0; i < poly.triangles.size(); i++) {
        out << 3 << " " << poly.triangles[i](0) << " " << poly.triangles[i](1) << " " << poly.triangles[i](2) << "\n";
    }

    for (int i = 0; i < poly.inner_verts.size(); i++) {
        out << poly.inner_verts[i](0) << " " << poly.inner_verts[i](1) << " " << 0.0 << "\n";
    }

    int bound_conds = 0;
    for (int i = 0; i < poly.triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {
            int v = poly.triangles[i](j);
            int w = (poly.triangles[i](j) + 1) % 3;
            if (poly.inner_verts[v](1) < -0.62 && poly.inner_verts[w](1) < -0.62) {
                bound_conds += 1;
            }

            if (dist(poly.inner_world[v], {0.3035, 0.4705}) < 0.05) {
                if (poly.bound_marks[v] && poly.bound_marks[w]) {
                    bound_conds += 1;
                }
            }
        }
    }

    out << bound_conds << " # Num. Boundary Conds.\n";
    for (int i = 0; i < poly.triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {
            int v = poly.triangles[i](j);
            int w = (poly.triangles[i](j) + 1) % 3;
            if (poly.inner_verts[v](1) < -0.62 && poly.inner_verts[w](1) < -0.62) {
                out << i << " " << j << " 40\n";
            }

            if (dist(poly.inner_world[v], {0.3035, 0.4705}) < 0.05) {
                if (poly.bound_marks[v] && poly.bound_marks[w]) {
                    out << i << " " << j << " 30\n";
                }
            }
        }
    }

    out.close();
}

// Begin the main program.
int main(int argc, char** argv) {
    RigidPoly obj = PolyReader::readSolid("/Users/teragion/Models/out.poly");
    Triangulate::triangulateSolid(obj);
    write_xda(obj);

    // Initialize libMesh and any dependent libraries
    LibMeshInit init(argc, argv);

    // This example requires a linear solver package.
    libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE, "--enable-petsc, --enable-trilinos, or --enable-eigen");

    // Initialize the cantilever mesh
    const unsigned int dim = 2;

    // Skip this 2D example if libMesh was compiled as 1D-only.
    libmesh_example_requires(dim <= LIBMESH_DIM, "2D support");

    // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
    libmesh_example_requires(false, "--enable-dirichlet");
#endif

    // Create a 2D mesh distributed across the default MPI communicator.
    Mesh mesh(init.comm(), dim);

    // 'MeshTool's::Generation::build_square(mesh, nx, ny, 0., 1., 0., 0.2, QUAD9);
    mesh.read("/Users/teragion/Models/out.xda");

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    // Declare the system and its variables.
    // Create a system named "Elasticity"
    // Declare the system and its variables.
    // Parse the input file
    GetPot infile("FEM.in");

    // Override input file arguments from the command line
    infile.parse_command_line(argc, argv);

    // Read in parameters from the input file
    const Real deltat = infile("deltat", 0.25);
    unsigned int n_timesteps = infile("n_timesteps", 1);

#ifdef LIBMESH_HAVE_EXODUS_API
    const unsigned int write_interval = infile("write_interval", 1);
#endif

    ElasticitySystem& system = equation_systems.add_system<ElasticitySystem>("Linear Elasticity");

    // Solve this as a time-dependent or steady system
    std::string time_solver = infile("time_solver", "newmark");

    ExplicitSystem* v_system;
    ExplicitSystem* a_system;

    if (time_solver == std::string("newmark")) {
        // Create ExplicitSystem to help output velocity
        v_system = &equation_systems.add_system<ExplicitSystem>("Velocity");
        v_system->add_variable("u_vel", FIRST, LAGRANGE);
        v_system->add_variable("v_vel", FIRST, LAGRANGE);
        v_system->add_variable("w_vel", FIRST, LAGRANGE);

        // Create ExplicitSystem to help output acceleration
        a_system = &equation_systems.add_system<ExplicitSystem>("Acceleration");
        a_system->add_variable("u_accel", FIRST, LAGRANGE);
        a_system->add_variable("v_accel", FIRST, LAGRANGE);
        a_system->add_variable("w_accel", FIRST, LAGRANGE);

        system.time_solver = std::make_unique<NewmarkSolver>(system);
    }

    else if (time_solver == std::string("euler")) {
        system.time_solver = std::make_unique<EulerSolver>(system);
        EulerSolver& euler_solver = cast_ref<EulerSolver&>(*(system.time_solver.get()));
        euler_solver.theta = infile("theta", 1.0);
    }

    else if (time_solver == std::string("euler2")) {
        system.time_solver = std::make_unique<Euler2Solver>(system);
        Euler2Solver& euler_solver = cast_ref<Euler2Solver&>(*(system.time_solver.get()));
        euler_solver.theta = infile("theta", 1.0);
    }

    else if (time_solver == std::string("steady")) {
        system.time_solver = std::make_unique<SteadySolver>(system);
        libmesh_assert_equal_to(n_timesteps, 1);
    } else
        libmesh_error_msg(std::string("ERROR: invalid time_solver ") + time_solver);

    ExplicitSystem& stress_system = equation_systems.add_system<ExplicitSystem>("StressSystem");
    unsigned int s00 = stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
    unsigned int s01 = stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
    unsigned int s11 = stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);

    unsigned int e00 = stress_system.add_variable("e_00", CONSTANT, MONOMIAL);
    unsigned int e01 = stress_system.add_variable("e_01", CONSTANT, MONOMIAL);
    unsigned int e11 = stress_system.add_variable("e_11", CONSTANT, MONOMIAL);

    // Initialize the system
    equation_systems.init();

    // Set the time stepping options
    system.deltat = deltat;

    // And the nonlinear solver options
    DiffSolver& solver = *(system.time_solver->diff_solver().get());
    solver.quiet = infile("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.max_nonlinear_iterations = infile("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance = infile("relative_step_tolerance", 1.e-3);
    solver.relative_residual_tolerance = infile("relative_residual_tolerance", 0.0);
    solver.absolute_residual_tolerance = infile("absolute_residual_tolerance", 0.0);

    // And the linear solver options
    solver.max_linear_iterations = infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance = infile("initial_linear_tolerance", 1.e-3);

    // Print information about the system to the screen.
    equation_systems.print_info();

    // If we're using EulerSolver or Euler2Solver, and we're using EigenSparseLinearSolver,
    // then we need to reset the EigenSparseLinearSolver to use SPARSELU because BICGSTAB
    // chokes on the Jacobian matrix we give it and Eigen's GMRES currently doesn't work.
    NewtonSolver* newton_solver = dynamic_cast<NewtonSolver*>(&solver);
    if (newton_solver && (time_solver == std::string("euler") || time_solver == std::string("euler2"))) {
#ifdef LIBMESH_HAVE_EIGEN_SPARSE
        LinearSolver<Number>& linear_solver = newton_solver->get_linear_solver();
        EigenSparseLinearSolver<Number>* eigen_linear_solver = dynamic_cast<EigenSparseLinearSolver<Number>*>(&linear_solver);

        if (eigen_linear_solver) eigen_linear_solver->set_solver_type(SPARSELU);
#endif
    }

    if (time_solver == std::string("newmark")) {
        NewmarkSolver* newmark = cast_ptr<NewmarkSolver*>(system.time_solver.get());
        newmark->compute_initial_accel();

        // Copy over initial velocity and acceleration for output.
        // Note we can do this because of the matching variables/FE spaces
        *(v_system->solution) = system.get_vector("_old_solution_rate");
        *(a_system->solution) = system.get_vector("_old_solution_accel");
    }

#ifdef LIBMESH_HAVE_EXODUS_API
    // Output initial state
    {
        std::ostringstream file_name;

        // We write the file in the ExodusII format.
        file_name << std::string("out.") + time_solver + std::string(".e-s.") << std::setw(3) << std::setfill('0') << std::right << 0;
        compute_stresses(equation_systems);

        ExodusII_IO(mesh).write_timestep(file_name.str(), equation_systems,
                                         1,  // This number indicates how many time steps
                                             // are being written to the file
                                         system.time);
    }
#endif  // #ifdef LIBMESH_HAVE_EXODUS_API

    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step = 0; t_step != n_timesteps; ++t_step) {
        // A pretty update message
        libMesh::out << "\n\nSolving time step " << t_step << ", time = " << system.time << std::endl;

        system.solve();
        compute_stresses(equation_systems);

        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();

        // Copy over updated velocity and acceleration for output.
        // Note we can do this because of the matching variables/FE spaces
        if (time_solver == std::string("newmark")) {
            *(v_system->solution) = system.get_vector("_old_solution_rate");
            *(a_system->solution) = system.get_vector("_old_solution_accel");
        }

#ifdef LIBMESH_HAVE_EXODUS_API
        // Write out this timestep if we're requested to
        if ((t_step + 1) % write_interval == 0) {
            std::ostringstream file_name;

            // We write the file in the ExodusII format.
            file_name << std::string("out.") + time_solver + std::string(".e-s.") << std::setw(3) << std::setfill('0') << std::right << t_step + 1;

            ExodusII_IO(mesh).write_timestep(file_name.str(), equation_systems,
                                             1,  // This number indicates how many time steps
                                                 // are being written to the file
                                             system.time);
        }
#endif  // #ifdef LIBMESH_HAVE_EXODUS_API
    }

    // All done.
    return 0;
}

void assemble_elasticity(EquationSystems& es, const std::string& libmesh_dbg_var(system_name)) {
    libmesh_assert_equal_to(system_name, "Elasticity");

    const MeshBase& mesh = es.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");

    const unsigned int u_var = system.variable_number("u");
    const unsigned int v_var = system.variable_number("v");

    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);

    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule(&qface);

    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kvu(Ke), Kvv(Ke);

    DenseSubVector<Number> Fu(Fe), Fv(Fe);

    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;

    SparseMatrix<Number>& matrix = system.get_system_matrix();

    for (const auto& elem : mesh.active_local_element_ptr_range()) {
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();

        fe->reinit(elem);

        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        Kuu.reposition(u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition(u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);

        Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);

        Fu.reposition(u_var * n_u_dofs, n_u_dofs);
        Fv.reposition(v_var * n_u_dofs, n_v_dofs);

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++) {
                    // Tensor indices
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 0, C_k = 0;

                    C_j = 0, C_l = 0;
                    Kuu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 0;
                    Kuu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 0, C_l = 1;
                    Kuu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 1;
                    Kuu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }

            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++) {
                    // Tensor indices
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 0, C_k = 1;

                    C_j = 0, C_l = 0;
                    Kuv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 0;
                    Kuv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 0, C_l = 1;
                    Kuv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 1;
                    Kuv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }

            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++) {
                    // Tensor indices
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 1, C_k = 0;

                    C_j = 0, C_l = 0;
                    Kvu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 0;
                    Kvu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 0, C_l = 1;
                    Kvu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 1;
                    Kvu(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }

            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++) {
                    // Tensor indices
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 1, C_k = 1;

                    C_j = 0, C_l = 0;
                    Kvv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 0;
                    Kvv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 0, C_l = 1;
                    Kvv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));

                    C_j = 1, C_l = 1;
                    Kvv(i, j) += JxW[qp] * (eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }
        }

        {
            for (auto side : elem->side_index_range())
                if (elem->neighbor_ptr(side) == nullptr) {
                    const std::vector<std::vector<Real>>& phi_face = fe_face->get_phi();
                    const std::vector<Real>& JxW_face = fe_face->get_JxW();

                    fe_face->reinit(elem, side);

                    if (mesh.get_boundary_info().has_boundary_id(elem, side, 2))  // Apply a traction on the right side
                    {
                        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                            for (unsigned int i = 0; i < n_v_dofs; i++) Fv(i) += JxW_face[qp] * (-1.) * phi_face[i][qp];
                        // for (unsigned int i = 0; i < n_u_dofs; i++) Fu(i) += JxW_face[qp] * (-0.1) * phi_face[i][qp];
                    }
                }
        }

        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        matrix.add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }
}

Real eval_elasticity_tensor(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    // Define the Poisson ratio
    const Real nu = 0.3;

    // Define Young's modulus
    const Real K = 800.0;

    // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
    const Real lambda_1 = K * nu / ((1. + nu) * (1. - 2. * nu));
    const Real lambda_2 = 0.5 * K / (1 + nu);

    // Define the Kronecker delta functions that we need here
    Real delta_ij = (i == j) ? 1. : 0.;
    Real delta_il = (i == l) ? 1. : 0.;
    Real delta_ik = (i == k) ? 1. : 0.;
    Real delta_jl = (j == l) ? 1. : 0.;
    Real delta_jk = (j == k) ? 1. : 0.;
    Real delta_kl = (k == l) ? 1. : 0.;

    return lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
}

/**
 * Compute the Cauchy stress for the current solution.
 */
void compute_stresses(EquationSystems& es) {
    // const Real young_modulus = es.parameters.get<Real>("young_modulus");
    // const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

    const MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // NonlinearImplicitSystem& system = es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Linear Elasticity");

    unsigned int displacement_vars[] = {system.variable_number("Ux"), system.variable_number("Uy")};
    const unsigned int u_var = system.variable_number("Ux");

    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&qrule);

    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

    // Also, get a reference to the ExplicitSystem
    ExplicitSystem& stress_system = es.get_system<ExplicitSystem>("StressSystem");
    const DofMap& stress_dof_map = stress_system.get_dof_map();
    unsigned int sigma_vars[] = {stress_system.variable_number("sigma_00"), stress_system.variable_number("sigma_01"), stress_system.variable_number("sigma_11")};

    unsigned int strain_vars[] = {stress_system.variable_number("e_00"), stress_system.variable_number("e_01"), stress_system.variable_number("e_11")};
    // Storage for the stress dof indices on each element
    std::vector<std::vector<dof_id_type>> dof_indices_var(system.n_vars());
    std::vector<dof_id_type> stress_dof_indices_var;
    std::vector<dof_id_type> strain_dof_indices_var;

    // To store the stress tensor on each element
    TensorValue<Number> elem_avg_stress_tensor;
    TensorValue<Number> elem_avg_strain_tensor;

    for (const auto& elem : mesh.active_local_element_ptr_range()) {
        for (unsigned int var = 0; var < 2; var++) dof_map.dof_indices(elem, dof_indices_var[var], displacement_vars[var]);

        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit(elem);

        // clear the stress tensor
        elem_avg_stress_tensor.zero();
        elem_avg_strain_tensor.zero();

        // printf("%lf\n", system.current_solution(dof_indices_var[0][0]));

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            TensorValue<Number> grad_u;
            // Row is variable u1, u2, or u3, column is x, y, or z
            for (unsigned int var_i = 0; var_i < 2; var_i++)
                for (unsigned int var_j = 0; var_j < 2; var_j++)
                    for (unsigned int j = 0; j < n_var_dofs; j++) grad_u(var_i, var_j) += dphi[j][qp](var_j) * system.current_solution(dof_indices_var[var_i][j]);

            TensorValue<Number> strain_tensor;
            for (unsigned int i = 0; i < 2; i++)
                for (unsigned int j = 0; j < 2; j++) {
                    strain_tensor(i, j) += 0.5 * (grad_u(i, j) + grad_u(j, i));

                    // effect of second order terms?
                    // for (unsigned int k = 0; k < 2; k++) {
                    //     strain_tensor(i, j) += 0.5 * grad_u(k, i) * grad_u(k, j);
                    //     if (grad_u(k, i) * grad_u(k, j) > 0.5) {
                    //         std::cout << "Alert!\n";
                    //     }
                    //     std::cout << grad_u(k, i) * grad_u(k, j) << " " << (grad_u(i, j) + grad_u(j, i)) << std::endl;
                    // }
                }

            // Define the deformation gradient
            auto F = grad_u;
            for (unsigned int var = 0; var < 2; var++) F(var, var) += 1.;

            TensorValue<Number> stress_tensor;
            for (unsigned int i = 0; i < 2; i++)
                for (unsigned int j = 0; j < 2; j++)
                    for (unsigned int k = 0; k < 2; k++)
                        for (unsigned int l = 0; l < 2; l++) stress_tensor(i, j) += eval_elasticity_tensor(i, j, k, l) * strain_tensor(k, l);

            // stress_tensor now holds the second Piola-Kirchoff stress (PK2) at point qp.
            // However, in this example we want to compute the Cauchy stress which is given by
            // 1/det(F) * F * PK2 * F^T, hence we now apply this transformation.
            auto Fdet = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
            stress_tensor = 1. / Fdet * F * stress_tensor * F.transpose();
            // stress_tensor = 1. / F.det() * F * stress_tensor * F.transpose();

            // We want to plot the average Cauchy stress on each element, hence
            // we integrate stress_tensor
            elem_avg_stress_tensor.add_scaled(stress_tensor, JxW[qp]);
            elem_avg_strain_tensor.add_scaled(strain_tensor, JxW[qp]);
        }

        // Get the average stress per element by dividing by volume
        elem_avg_stress_tensor /= elem->volume();

        // load elem_sigma data into stress_system
        unsigned int var_index = 0;
        for (unsigned int i = 0; i < 2; i++)
            for (unsigned int j = i; j < 2; j++) {
                stress_dof_map.dof_indices(elem, stress_dof_indices_var, sigma_vars[var_index]);
                stress_dof_map.dof_indices(elem, strain_dof_indices_var, strain_vars[var_index]);

                // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
                // one dof index per variable
                dof_id_type dof_index = stress_dof_indices_var[0];

                if ((stress_system.solution->first_local_index() <= dof_index) && (dof_index < stress_system.solution->last_local_index())) {
                    stress_system.solution->set(dof_index, elem_avg_stress_tensor(i, j));
                }

                dof_index = strain_dof_indices_var[0];

                if ((stress_system.solution->first_local_index() <= dof_index) && (dof_index < stress_system.solution->last_local_index())) {
                    stress_system.solution->set(dof_index, elem_avg_strain_tensor(i, j));
                }

                var_index++;
            }
    }

    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
}
