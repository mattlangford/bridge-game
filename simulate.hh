#pragma once

// Seems like a bug in Eigen requires us to ignore this warning
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-anon-enum-enum-conversion"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

#include "mesh.hh"

#include <GLFW/glfw3.h>

using DMatrix = Eigen::Matrix3d;
using BMatrix = Eigen::Matrix<double, 3, 6>;

static constexpr double kDt = 1.0 / 500.0;

//
// #############################################################################
//

DMatrix generate_D()
{
    // As a proof of concept I'm just going to hardcode these
    const double E = 3.7 * 1E5; // youngs modulus N/m^2 (for brick) (slightly adjusted)
    const double v = 0.1; // poissons ratio (also for brick)

    // Comes from [1] 4.14
    Eigen::Matrix3d d;
    // clang-format off
    d << 1.0,   v, 0.0,
           v, 1.0, 0.0,
         0.0, 0.0, 0.5f * (1.0 - v);
    // clang-format on

    return d * E / (1 - (v * v));
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
    void step(const double dt)
    {
        const auto start = std::chrono::high_resolution_clock::now();
        const size_t num_steps = dt / kDt;
        for (size_t i = 0; i < num_steps; ++i) {
            step();
        }

        std::cout << "Performed " << num_steps << " steps in "
                  << std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::high_resolution_clock::now() - start)
                         .count()
                  << "us\n";
    }

    void step()
    {
        const size_t vertex_count = U_.size();

        // From [2] Table 9.3 A.4, we'll define some constants to help out later
        constexpr double kAlpha = 0.5f;
        constexpr double kBeta = 0.25f * (0.5f + kAlpha) * (0.5f + kAlpha);
        constexpr double kA0 = 1.0f / (kBeta * kDt * kDt);
        constexpr double kA1 = kAlpha / (kBeta * kDt);
        constexpr double kA2 = 1.0f / (kBeta * kDt);
        constexpr double kA3 = 1.0f / (2.0f * kBeta) - 1.0f;
        constexpr double kA4 = kAlpha / kBeta - 1.0f;
        constexpr double kA5 = (kDt / 2.0f) * (kAlpha / kBeta - 2.0f);
        constexpr double kA6 = kDt * (1.0f - kAlpha);
        constexpr double kA7 = kDt * kAlpha;

        // Generate the mass matrix
        if (!mass_) {
            Eigen::MatrixXd& M = mass_.emplace(vertex_count, vertex_count);
            M.setZero();
            for (size_t i = 0; i < mesh_.mass.size(); ++i) {
                const size_t index = vertex_to_u_[i];
                if (index < U_.size()) {
                    M(index, index) = mesh_.mass[i];
                }
            }
        }
        Eigen::MatrixXd& M = *mass_;

        // Generate arbitrary dampening matrix
        if (!damping_) {
            // Not really sure what this should be
            constexpr double kDampingFactor = 10.0;
            Eigen::MatrixXd& C = damping_.emplace(vertex_count, vertex_count);
            C.setZero();
            C.diagonal().fill(kDampingFactor);
        }
        Eigen::MatrixXd& C = *damping_;

        if (!K_solver_) {
            // Generate the global stiffness matrix
            Eigen::SparseMatrix<double> K(vertex_count, vertex_count);
            K.setZero();
            K.reserve(6 * 6 * vertex_count);
            for (size_t i = 0; i < mesh_.triangles.size(); ++i) {
                // No processing needed if the triangle is fixed
                if (fixed_triangles_[i])
                    continue;

                auto& triangle = mesh_.triangles[i];

                static constexpr double kThickness = 1; // meters

                const BMatrix B = generate_B(triangle);
                Eigen::Matrix<double, 6, 6> k = kThickness * area(triangle) * B.transpose() * generate_D() * B;

                for (const auto& [local, global] : generate_local_to_global_mapping(triangle)) {
                    const auto& [local_row, local_col] = local;
                    const auto& [global_row, global_col] = global;

                    // We get the global indices above, but since the K matrix only includes dynamic vertices we'll
                    // need to convert one more time. If the global row/col map to a fixed vertex, we'll ignore it
                    const size_t dynamic_row = vertex_to_u_[global_row];
                    const size_t dynamic_col = vertex_to_u_[global_col];
                    if (dynamic_row >= U_.size() || dynamic_col >= U_.size()) {
                        continue;
                    }

                    K.coeffRef(dynamic_row, dynamic_col) += k(local_row, local_col);
                }
            }

            // [2] Table 9.3 A.5
            K = K + kA0 * M + kA1 * C;

            // [2] Table 9.3 A.6
            K_solver_.emplace();
            K_solver_->compute(K);
            if (K_solver_->info() != Eigen::Success) {
                throw std::runtime_error("Decomposition Failed!");
            }
        }

        Eigen::VectorXd gravity = Eigen::VectorXd::Zero(vertex_count);
        for (size_t i = 1; i < vertex_count; i += 2) {
            gravity[i] = -9.8 * M(i, i);
        }

        // [2] Table 9.3 B.1
        Eigen::VectorXd R_hat = gravity
            + M * (kA0 * U_ + kA2 * U_vel_ + kA3 * U_accel_)
            + C * (kA1 * U_ + kA4 * U_vel_ + kA5 * U_accel_);

        // [2] Table 9.3 B.2
        Eigen::VectorXd U_next = K_solver_->solve(R_hat);
        if (K_solver_->info() != Eigen::Success) {
            throw std::runtime_error("Solving Failed!");
        }

        // [2] Table 9.3 B.3
        Eigen::VectorXd U_accel_next = kA0 * (U_next - U_) - kA2 * U_vel_ - kA3 * U_accel_;
        Eigen::VectorXd U_vel_next = U_vel_ + kA6 * U_accel_ + kA7 * U_accel_next;

        // Now we can update our internal state!
        U_ = std::move(U_next);
        U_vel_ = std::move(U_vel_next);
        U_accel_ = std::move(U_accel_next);

        // To avoid numerical issues with things moving too quickly, we'll impose a max velocity
        for (size_t i = 0; i < U_vel_.size(); ++i) {
            double& velocity = U_vel_[i];

            constexpr double kTerminalVelocity = 30.0;
            velocity = std::clamp(velocity, -kTerminalVelocity, kTerminalVelocity);
        }
    }

    double compute_stress(const Triangle& triangle) const
    {
        Eigen::Matrix<double, 6, 1> displacements;
        for (size_t i = 0; i < triangle.indices.size(); ++i) {
            displacements[2 * i] = get_displacement(triangle.indices[i]);
            displacements[2 * i + 1] = get_displacement(triangle.indices[i] + 1);
        }

        Eigen::Vector3d stresses = generate_D() * generate_B(triangle) * displacements;
        return stresses.norm();
    }

    void draw() const
    {
        glBegin(GL_TRIANGLES);

        for (size_t i = 0; i < mesh_.triangles.size(); ++i) {
            const Triangle& triangle = mesh_.triangles[i];

            if (fixed_triangles_[i]) {
                glColor3f(0.f, 0.f, 1.f);
            } else {
                constexpr double kMaxStress = 10000;
                const double stress = compute_stress(triangle);

                float red = static_cast<float>(stress / kMaxStress); // 0 when stress is 0, 1 when stress is high
                float green = static_cast<float>((kMaxStress - stress) / kMaxStress); // 1 when stress is 0, 0 when stress is high
                glColor3f(std::clamp(red, 0.f, 1.f), std::clamp(green, 0.f, 1.f), 0.0f);
            }

            constexpr size_t kPixelSize = 10;
            const auto x = [&](uint8_t i) { return kPixelSize * get_coordinate(triangle.indices[i]); };
            const auto y = [&](uint8_t i) { return kPixelSize * get_coordinate(triangle.indices[i] + 1); };

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

        U_ = Eigen::VectorXd::Zero(num_dynamic);
        U_vel_ = Eigen::VectorXd::Zero(num_dynamic);
        U_accel_ = Eigen::VectorXd::Zero(num_dynamic);

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

        K_solver_ = std::nullopt;
        mass_ = std::nullopt;
        damping_ = std::nullopt;
    }

private:
    double get_displacement(size_t index) const
    {
        const size_t u_index = vertex_to_u_[index];
        return u_index >= U_.size() ? 0.0f : U_[u_index];
    }
    double get_coordinate(size_t index) const
    {
        return mesh_.vertices[index] + get_displacement(index);
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

        const double det_j = x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3);

        BMatrix b;
        // clang-format off
        b << y(2, 3),     0.0, y(3, 1),     0.0, y(1, 2),     0.0,
                 0.0, x(3, 2),     0.0, x(1, 3),     0.0, x(2, 1),
             x(3, 2), y(2, 3), x(1, 3), y(3, 1), x(2, 1), y(1, 2);
        // clang-format on
        return b / det_j;
    }

    double area(const Triangle& triangle) const
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
    Eigen::VectorXd U_;
    Eigen::VectorXd U_vel_;
    Eigen::VectorXd U_accel_;

    // Mesh we've been blessed with
    Mesh mesh_;

    std::optional<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> K_solver_;
    std::optional<Eigen::MatrixXd> mass_;
    std::optional<Eigen::MatrixXd> damping_;

    // Since we'll only generate displacements for non-fixed vertices, we'll need to store a mapping between mesh
    // vertices and displacement/velocity/accel vectors
    std::vector<size_t> vertex_to_u_;

    // Generated so we know which triangles are completely fixed
    std::vector<bool> fixed_triangles_;
};

#pragma clang diagnostic pop
