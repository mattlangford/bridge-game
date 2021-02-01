#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>

#include "common/mesh.hh"

using DMatrix = Eigen::Matrix3d;
using BMatrix = Eigen::Matrix<double, 3, 6>;
using GlobalKMatrix = Eigen::SparseMatrix<double>;
using LocalKMatrix = Eigen::Matrix<double, 6, 6>;

//
// #############################################################################
//

DMatrix generate_D();

//
// #############################################################################
//

Eigen::MatrixXd generate_mass_matrix(size_t num_displacements, const std::vector<size_t> &vertex_to_displacements,
                                     const std::vector<double> &mass);

//
// #############################################################################
//

Eigen::MatrixXd generate_damping_matrix(size_t vertex_count);

//
// #############################################################################
//

struct State {
    // Displacements of the dynamic coordinates (along with per-node velocity and
    // accel)
    Eigen::VectorXd displacements;
    Eigen::VectorXd velocities;
    Eigen::VectorXd accelerations;
};

struct Cache {
    Eigen::MatrixXd mass;
    Eigen::MatrixXd damping;
    Eigen::SimplicialLDLT<GlobalKMatrix> K_solver;

    // Since we'll only generate displacements for non-fixed vertices, we'll need
    // to store a mapping between mesh vertices and displacement/velocity/accel
    // vectors
    std::vector<size_t> vertex_to_displacements;

    // Generated so we know which triangles are completely fixed
    std::vector<bool> fixed_triangles;
};

struct SimulationContext {
    SimulationContext(common::Mesh mesh_);

    State state;
    common::Mesh mesh;
    Cache cache;

    inline double get_displacement(size_t index) const {
        const size_t u_index = cache.vertex_to_displacements[index];
        return u_index >= state.displacements.size() ? 0.0f : state.displacements[u_index];
    }

    inline double get_coordinate(size_t index) const { return mesh.vertices[index] + get_displacement(index); }
};

struct MeshStepper {
    // From [2] Table 9.3 A.4, we'll define some constants to help out later
    static constexpr double kDt = 1.0 / 250.0;
    static constexpr double kAlpha = 0.5f;
    static constexpr double kBeta = 0.25f * (0.5f + kAlpha) * (0.5f + kAlpha);
    static constexpr double kA0 = 1.0f / (kBeta * kDt * kDt);
    static constexpr double kA1 = kAlpha / (kBeta * kDt);
    static constexpr double kA2 = 1.0f / (kBeta * kDt);
    static constexpr double kA3 = 1.0f / (2.0f * kBeta) - 1.0f;
    static constexpr double kA4 = kAlpha / kBeta - 1.0f;
    static constexpr double kA5 = (kDt / 2.0f) * (kAlpha / kBeta - 2.0f);
    static constexpr double kA6 = kDt * (1.0f - kAlpha);
    static constexpr double kA7 = kDt * kAlpha;

    static State step(SimulationContext &context);
    std::vector<double> evaluate_stresses() const;
};

struct TriangleStressHelpers {
    const common::Triangle &triangle;
    const SimulationContext &context;

    // These both assume the same formatting as found in [1] 4.71 and 4.72. They
    // are 1 indexed, and assume indices are ordered x0, y0, x1, y1, ...
    double x(size_t i, size_t j) const;
    double y(size_t i, size_t j) const;

    struct LocalToGlobalKMatrixMapping {
        std::pair<size_t, size_t> local_row_col;
        std::pair<size_t, size_t> global_row_col;
    };

    std::vector<LocalToGlobalKMatrixMapping> generate_local_to_global_mapping() const;

    void populate_local_stiffness(GlobalKMatrix &global_k) const;

    BMatrix generate_B() const;

    double area() const;
};

GlobalKMatrix generate_global_stiffness_matrix(const SimulationContext &context);

//
// #############################################################################
//

class Simulator {
   public:
    void step(const double dt);

    void draw() const;

    void set_mesh(common::Mesh mesh);

   private:
    void destroy_stressful_triangles();

    Eigen::Vector3d compute_stress(const common::Triangle &triangle) const;

   private:
    std::unique_ptr<SimulationContext> context;
};
