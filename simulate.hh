#pragma once

// Seems like a bug in Eigen requires us to ignore this warning
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-anon-enum-enum-conversion"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <iostream>
#include <optional>
#include <vector>

#include "mesh.hh"

#include <GLFW/glfw3.h>

using DMatrix = Eigen::Matrix3f;
using BMatrix = Eigen::Matrix<float, 3, 6>;

static constexpr float kDt = 1.f / 60.f;

//
// #############################################################################
//

DMatrix generate_D()
{
    // As a proof of concept I'm just going to hardcode these
    const float E = 3.7 * 1E10; // youngs modulus N/m^2 (for brick)
    const float v = 0.1; // poissons ratio (also for brick)

    // Comes from [1] 4.14
    Eigen::Matrix3f d;
    // clang-format off
    d << 1.f,   v, 0.f,
           v, 1.f, 0.f,
         0.f, 0.f, 0.5f * (1.f - v);
    // clang-format on

    return E / (1.f - (v * v)) * d;
}

//
// #############################################################################
//

struct LocalToGlobalMapping {
    std::pair<size_t, size_t> local_row_col;
    std::pair<size_t, size_t> global_row_col;
};

std::vector<LocalToGlobalMapping> generate_local_to_global_mapping(const Triangle& triangle)
{
    std::vector<LocalToGlobalMapping> mapping;
    mapping.reserve(6 * 6);

    for (size_t row = 0; row < 3; ++row) {
        for (size_t col = 0; col < 3; ++col) {
            // Map the X values
            auto& x_map = mapping.emplace_back();
            x_map.local_row_col = { 2 * row, 2 * col };
            x_map.global_row_col = { triangle.indices[row], triangle.indices[col] };

            // Map the Y values
            auto& y_map = mapping.emplace_back();
            y_map.local_row_col = { 2 * row + 1, 2 * col + 1 };
            y_map.global_row_col = { triangle.indices[row] + 1, triangle.indices[col] + 1 };
        }
    }

    return mapping;
}

//
// #############################################################################
//

class Simlator {
public:
    void step()
    {
        const size_t vertex_count = U_.size();

        // From [2] Table 9.3 A.4, we'll define some constants to help out later
        constexpr float kAlpha = 0.5f;
        constexpr float kBeta = 0.25f * (0.5f + kAlpha) * (0.5f + kAlpha);
        constexpr float kA0 = 1.0f / (kBeta * kDt * kDt);
        constexpr float kA1 = kAlpha / (kBeta * kDt);
        constexpr float kA2 = 1.0f / (kBeta * kDt);
        constexpr float kA3 = 1.0f / (2.0f * kBeta) - 1.0f;
        constexpr float kA4 = kAlpha / kBeta - 1.0f;
        constexpr float kA5 = (kDt / 2.0f) * (kAlpha / kBeta - 2.0f);
        constexpr float kA6 = kDt * (1.0f - kAlpha);
        constexpr float kA7 = kDt * kAlpha;

        Eigen::SparseMatrix<float> K(vertex_count, vertex_count);
        K.setZero();
        K.reserve(6 * 6 * vertex_count);

        Eigen::MatrixXf M(vertex_count, vertex_count);
        M.setZero();
        for (size_t i = 0; i < mesh_.mass.size(); ++i) {
            const size_t index = vertex_to_u_[i];
            if (index < U_.size()) {
                M(index, index) = mesh_.mass[i];
            }
        }

        for (size_t i = 0; i < mesh_.triangles.size(); ++i) {
            // No processing needed if the triangle is fixed
            if (fixed_triangles_[i])
                continue;

            auto& triangle = mesh_.triangles[i];

            static constexpr float kThickness = 1; // meters

            const BMatrix B = generate_B(triangle);
            Eigen::Matrix<float, 6, 6> k = kThickness * area(triangle) * B.transpose() * generate_D() * B;

            for (const auto& [local, global] : generate_local_to_global_mapping(triangle)) {
                const auto& [local_row, local_col] = local;
                const auto& [global_row, global_col] = global;

                // We get the global indices above, but since the K matrix only includes dynamic vertices we'll need to
                // convert one more time. If the global row/col map to a fixed vertex, we'll ignore this entry
                const size_t dynamic_row = vertex_to_u_[global_row];
                const size_t dynamic_col = vertex_to_u_[global_col];
                if (dynamic_row >= U_.size() || dynamic_col >= U_.size()) {
                    continue;
                }

                K.coeffRef(dynamic_row, dynamic_col) += k(local_row, local_col);
            }
        }

        // [2] Table 9.3 A.5 (with C = 0)
        K = K + kA0 * M;

        // [2] Table 9.3 A.6
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
        solver.compute(K);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Decomposition Failed!");
        }

        // [2] Table 9.3 B.1 (again with C = 0)
        Eigen::VectorXf R_hat = M * (kA0 * U_ + kA2 * U_vel_ + kA3 * U_accel_);

        // [2] Table 9.3 B.2
        Eigen::VectorXf U_next = solver.solve(R_hat);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving Failed!");
        }

        // [2] Table 9.3 B.3
        Eigen::VectorXf U_accel_next = kA0 * (U_next - U_) - kA2 * U_vel_ - kA3 * U_accel_;
        Eigen::VectorXf U_vel_next = U_vel_ + kA6 * U_accel_ + kA7 * U_accel_next;

        // Now we can update our internal state!
        U_ = std::move(U_next);
        U_vel_ = std::move(U_vel_next);
        U_accel_ = std::move(U_accel_next);
    }

    void draw() const
    {
        glBegin(GL_TRIANGLES);

        Metadata last_metadata;
        for (const Triangle& triangle : mesh_.triangles) {
            const auto x = [&](uint8_t i) { return get_coordinate(triangle.indices[i]); };
            const auto y = [&](uint8_t i) { return get_coordinate(triangle.indices[i] + 1); };

            glVertex2f(x(0), y(0));
            glVertex2f(x(1), y(1));
            glVertex2f(x(2), y(2));
        }

        glEnd();
    }

    void set_mesh(Mesh mesh)
    {
        mesh_ = std::move(mesh);

        const size_t vertex_count = mesh_.vertices.size();

        size_t num_dynamic = 0;
        vertex_to_u_.resize(vertex_count, -1);
        for (size_t i = 0; i < vertex_count; ++i) {
            if (mesh_.fixed[i])
                continue;

            vertex_to_u_[i] = num_dynamic; // store before incrementing
            num_dynamic++;
        }

        U_ = Eigen::VectorXf::Zero(num_dynamic);
        U_vel_ = Eigen::VectorXf::Zero(num_dynamic);
        U_accel_ = Eigen::VectorXf::Zero(num_dynamic);

        // Apply gravity (only to the Y coordinates)
        for (size_t i = 1; i < U_accel_.size(); i += 2) {
            U_accel_[i] = -9.8;
        }

        fixed_triangles_.reserve(vertex_count);
        for (const auto& triangle : mesh_.triangles) {
            // Fixed until proven otherwise
            bool fixed = true;

            for (size_t index : triangle.indices) {
                if (!mesh_.fixed[index]) {
                    // At least one index isn't fixed!
                    fixed = false;
                    break;
                }
            }

            fixed_triangles_.push_back(fixed);
        }
    }

private:
    float get_coordinate(size_t index) const
    {
        const size_t u_index = vertex_to_u_[index];
        const float displacement = u_index >= U_.size() ? 0.0f : U_[u_index];
        return mesh_.vertices[index] + displacement;
    }

    BMatrix generate_B(const Triangle& triangle) const
    {
        auto x = [&](size_t i, size_t j) {
            // we need to subtract one to get the indexing to match [1]
            return get_coordinate(triangle.indices[i - 1]) - get_coordinate(triangle.indices[j - 1]);
        };
        auto y = [&](size_t i, size_t j) {
            // we need to subtract one to get the indexing to match [1]
            // we also have to add one since the y indices come after the x ones
            return get_coordinate(triangle.indices[i - 1] + 1) - get_coordinate(triangle.indices[j - 1] + 1);
        };

        const float det_j = x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3);

        BMatrix b;
        // clang-format off
        b << y(2, 3),     0.f, y(3, 1),     0.f, y(1, 2),     0.f,
                 0.f, x(3, 2),     0.f, x(1, 3),     0.f, x(2, 1),
             x(3, 2), y(2, 3), x(1, 3), y(3, 1), x(2, 1), y(1, 2);
        // clang-format on
        return b / det_j;
    }

    float area(const Triangle& triangle) const
    {
        auto x = [&](size_t i, size_t j) {
            // we need to subtract one to get the indexing to match [1]
            return get_coordinate(triangle.indices[i - 1]) - get_coordinate(triangle.indices[j - 1]);
        };
        auto y = [&](size_t i, size_t j) {
            // we need to subtract one to get the indexing to match [1]
            // we also have to add one since the y indices come after the x ones
            return get_coordinate(triangle.indices[i - 1] + 1) - get_coordinate(triangle.indices[j - 1] + 1);
        };

        // From [1] 4.70, I'm keeping the indices one-indexed like they do
        return 0.5f * abs(x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3));
    }

private:
    // Displacements of the dynamic coordinates (along per-node velocity and accel)
    Eigen::VectorXf U_;
    Eigen::VectorXf U_vel_;
    Eigen::VectorXf U_accel_;

    // Mesh we've been blessed with
    Mesh mesh_;

    // Since we'll only generate displacements for non-fixed vertices, we'll need to store a mapping between mesh
    // vertices and displacement/velocity/accel vectors
    std::vector<size_t> vertex_to_u_;

    // Generated so we know which triangles are completely fixed
    std::vector<bool> fixed_triangles_;
};

#pragma clang diagnostic pop
