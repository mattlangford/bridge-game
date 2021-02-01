#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "common/mesh.hh"

using DMatrix = Eigen::Matrix3d;
using BMatrix = Eigen::Matrix<double, 3, 6>;
using GlobalKMatrix = Eigen::SparseMatrix<double>;
using LocalKMatrix = Eigen::Matrix<double, 6, 6>;

///
/// @brief Current state of the system. This will hold per-point displacements, velocities, and accelerations
///
struct State {
    Eigen::VectorXd displacements;
    Eigen::VectorXd velocities;
    Eigen::VectorXd accelerations;
};

///
/// @brief Precomputed cache to speed up the main dynamics loop
///
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

///
/// @brief Full context used for simulation
///
struct SimulationContext {
    SimulationContext(common::Mesh mesh_);

    common::Mesh mesh;
    State state;
    Cache cache;

    inline double get_displacement(size_t index) const {
        const size_t u_index = cache.vertex_to_displacements[index];
        return u_index >= static_cast<size_t>(state.displacements.size()) ? 0.0f : state.displacements[u_index];
    }

    inline double get_coordinate(size_t index) const { return mesh.vertices[index] + get_displacement(index); }
};
