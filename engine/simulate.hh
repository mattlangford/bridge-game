#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

#include "common/mesh.hh"
#include "engine/context.hh"

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

struct MeshStepper {
    // From [2] Table 9.3 A.4, we'll define some constants to help out later
    static constexpr double kDt = 1.0 / 500.0;
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

    void set_mesh(common::Mesh mesh);
    const SimulationContext *simulation_context() const;

   private:
    void destroy_stressful_triangles();

    Eigen::Vector3d compute_stress(const common::Triangle &triangle) const;

   private:
    std::unique_ptr<SimulationContext> context;
};
