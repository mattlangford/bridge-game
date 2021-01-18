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

using DMatrix = Eigen::Matrix3f;
using BMatrix = Eigen::Matrix<float, 3, 6>;

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

    return (E / 1.f - (v * v)) * d;
}

//
// #############################################################################
//

BMatrix generate_B(const Triangle& triangle, const Mesh& mesh)
{
    auto x = [&](size_t i, size_t j) {
        // to get the indexing to match [1]
        i -= 1;
        j -= 1;
        return mesh.vertices[triangle.indices[i]] - mesh.vertices[triangle.indices[j]];
    };
    auto y = [&](size_t i, size_t j) {
        // to get the indexing to match [1]
        i -= 1;
        j -= 1;
        return mesh.vertices[triangle.indices[i] + 1] - mesh.vertices[triangle.indices[j] + 1];
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
        const size_t vertex_count = mesh_.vertices.size();

        // From [2] Table 9.1 A.3, we'll define some constants to help out later
        constexpr float kDt = 1.f / 60.f;
        constexpr float kA0 = 1.f / (kDt * kDt);
        constexpr float kA1 = 1.f / (2.f * kDt);
        constexpr float kA2 = 2.f * kA0;
        constexpr float kA3 = 1.f / kA2;

        // Out of laze, I'm going to store the full K matrix!
        Eigen::SparseMatrix<float> K(vertex_count, vertex_count);
        K.setZero();
        K.reserve(6 * 6 * vertex_count);

        // Mass of each node, for now it'll be just a third of each triangles mass
        Eigen::MatrixXf M(vertex_count, vertex_count);
        M.setZero();

        for (Triangle& triangle : mesh_.triangles) {
            if (triangle.metadata.fixed) {
                continue;
            }

            static constexpr float kThickness = 1; // meters

            const BMatrix B = generate_B(triangle, mesh_);
            Eigen::Matrix<float, 6, 6> k = kThickness * area(triangle, mesh_) * B.transpose() * generate_D() * B;

            for (const auto& [local, global] : generate_local_to_global_mapping(triangle)) {
                const auto& [local_row, local_col] = local;
                const auto& [global_row, global_col] = global;

                K.coeffRef(global_row, global_col) += k(local_row, local_col);
            }

            const float triangle_mass_sixth = triangle.metadata.mass / 6.f;
            for (size_t index : triangle.indices) {
                M.diagonal()[index] += triangle_mass_sixth;
                M.diagonal()[index + 1] += triangle_mass_sixth;
            }
        }

        // Now that we have effective masses for each node, we can calculate nodal gravity forces
        Eigen::VectorXf gravity(vertex_count);
        gravity.setZero();
        for (size_t node = 1; node < vertex_count; node += 2) // only iterate over the Y values
        {
            constexpr float kA = -9.8; // acceleration due to gravity (m/s2)
            gravity[node] = M.diagonal()[node] * kA;
        }

        Eigen::VectorXf U = Eigen::Map<Eigen::VectorXf>(mesh_.vertices.data(), mesh_.vertices.size());

        // Eq from [2] Table 9.1 B.1
        Eigen::VectorXf R_hat = gravity;
        R_hat -= K * U;
        R_hat -= (kA2 * M) * U;
        R_hat -= (kA0 * M) * U_prev_;
        for (size_t i = 0; i < vertex_count; ++i) {
            const float mass = M.diagonal()[i];
            if (mass <= 0.f) {
                continue;
            }

            std::cout << R_hat[i] << "\n";
            mesh_.vertices[i] = R_hat[i] / mass;
        }
        exit(1);
    }

    // void step()
    // {
    //     if (vertex_velocities_.size() != mesh_.vertices.size()) {
    //         throw std::runtime_error("Invalid mesh!");
    //     }

    //     std::vector<bool> fixed_mask(mesh_.vertices.size(), false);
    //     for (Triangle& triangle : mesh_.triangles) {
    //         if (!triangle.metadata.fixed) {
    //             continue;
    //         }

    //         for (size_t index : triangle.indices) {
    //             fixed_mask[index] = true;
    //             fixed_mask[index + 1] = true;
    //         }
    //     }

    //     for (size_t i = 0; i < mesh_.vertices.size(); ++i) {
    //         if (i % 2 == 0)
    //             continue;
    //         if (fixed_mask[i])
    //             continue;

    //         vertex_velocities_[i] -= 0.1;
    //         mesh_.vertices[i] += vertex_velocities_[i];
    //     }
    // }

    void draw() const
    {
        draw_mesh(mesh_);
    }

    void set_mesh(Mesh mesh)
    {
        mesh_ = std::move(mesh);

        U_prev_ = Eigen::VectorXf::Zero(mesh_.vertices.size());
        for (size_t i = 0; i < mesh_.vertices.size(); ++i) {
            U_prev_[i] = mesh_.vertices[i];
        }
    }

private:
    Eigen::VectorXf U_prev_;
    Mesh mesh_;
};

#pragma clang diagnostic pop
