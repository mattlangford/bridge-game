#include "engine/simulate.hh"

#include <GLFW/glfw3.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

#include "common/config.hh"

static Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

//
// #############################################################################
//

DMatrix generate_D() {
    // As a proof of concept I'm just going to hardcode these
    const double E = 30 * 1E9;  // youngs modulus N/m^2 (for concrete) (slightly adjusted)
    const double v = 0.25;      // poissons ratio (also for concrete)

    // Comes from [1] 4.14
    Eigen::Matrix3d d;
    // clang-format off
    d << 1.0,   v, 0.0,
           v, 1.0, 0.0,
         0.0, 0.0, (1.0 - v) / 2.0;
    // clang-format on

    return d * E / (1.0 - v * v);
}

//
// #############################################################################
//

Eigen::MatrixXd generate_mass_matrix(size_t num_displacements, const std::vector<size_t> &vertex_to_displacements,
                                     const std::vector<double> &mass) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(num_displacements, num_displacements);
    M.setZero();
    for (size_t i = 0; i < mass.size(); ++i) {
        const size_t index = vertex_to_displacements[i];
        if (index < num_displacements) {
            M(index, index) += mass[i];
        }
    }
    return M;
}

//
// #############################################################################
//

Eigen::MatrixXd generate_damping_matrix(size_t vertex_count) {
    // Not really sure what this should be
    constexpr double kDampingFactor = 0.0;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(vertex_count, vertex_count);
    C.diagonal().fill(kDampingFactor);
    return C;
}

//
// #############################################################################
//

double TriangleStressHelpers::x(size_t i, size_t j) const {
    auto i_coord = context.get_coordinate(triangle.indices[i - 1]);
    auto j_coord = context.get_coordinate(triangle.indices[j - 1]);
    return i_coord - j_coord;
}
double TriangleStressHelpers::y(size_t i, size_t j) const {
    auto i_coord = context.get_coordinate(triangle.indices[i - 1] + 1);
    auto j_coord = context.get_coordinate(triangle.indices[j - 1] + 1);
    return i_coord - j_coord;
}

//
// #############################################################################
//

auto TriangleStressHelpers::generate_local_to_global_mapping() const -> std::vector<LocalToGlobalKMatrixMapping> {
    std::vector<LocalToGlobalKMatrixMapping> mapping;
    mapping.reserve(6 * 6);

    for (size_t row = 0; row < 3; ++row) {
        for (size_t col = 0; col < 3; ++col) {
            // Map the X values
            auto &x_map = mapping.emplace_back();
            x_map.local_row_col = {2 * row, 2 * col};
            x_map.global_row_col = {triangle.indices[row], triangle.indices[col]};

            // Map the Y values
            auto &y_map = mapping.emplace_back();
            y_map.local_row_col = {2 * row + 1, 2 * col + 1};
            y_map.global_row_col = {triangle.indices[row] + 1, triangle.indices[col] + 1};
        }
    }

    return mapping;
}

//
// #############################################################################
//

void TriangleStressHelpers::populate_local_stiffness(GlobalKMatrix &global_k) const {
    static constexpr double kThickness = common::kBlockSize;  // meters

    const BMatrix B = generate_B();
    const size_t num_displacements = context.state.displacements.size();

    // [1] 4.66
    LocalKMatrix local_k = kThickness * area() * B.transpose() * generate_D() * B;
    // for (const auto& [local, global] : generate_local_to_global_mapping()) {
    //     const auto& [local_row, local_col] = local;
    //     const auto& [global_row, global_col] = global;

    //     // We get the global indices above, but since the K matrix only
    //     includes dynamic vertices we'll
    //     // need to convert one more time. If the global row/col map to a
    //     fixed vertex, we'll ignore it const auto& vertex_to_displacements =
    //     context.cache.vertex_to_displacements; const size_t dynamic_row =
    //     vertex_to_displacements[global_row]; const size_t dynamic_col =
    //     vertex_to_displacements[global_col]; std::cout << "local: " <<
    //     local_row << ", " << local_col << ", "; std::cout << "global: " <<
    //     global_row << ", " << global_col << ", "; std::cout << "dynamic: " <<
    //     dynamic_row << ", " << dynamic_col << ", "; std::cout << "\n"; if
    //     (dynamic_row < num_displacements && dynamic_col < num_displacements)
    //     {
    //         global_k.coeffRef(dynamic_row, dynamic_col) += local_k(local_row,
    //         local_col);
    //     }
    // }
    auto to_global = [this](size_t i) {
        if (i % 2 == 0) {
            return triangle.indices[i / 2];
        } else {
            return triangle.indices[i / 2] + 1;
        }
    };
    for (size_t local_row = 0; local_row < 6; ++local_row) {
        for (size_t local_col = 0; local_col < 6; ++local_col) {
            const size_t global_row = to_global(local_row);
            const size_t global_col = to_global(local_col);

            const auto &vertex_to_displacements = context.cache.vertex_to_displacements;
            const size_t dynamic_row = vertex_to_displacements[global_row];
            const size_t dynamic_col = vertex_to_displacements[global_col];
            if (dynamic_row < num_displacements && dynamic_col < num_displacements) {
                global_k.coeffRef(dynamic_row, dynamic_col) += local_k(local_row, local_col);
            }
        }
    }
}

//
// #############################################################################
//

BMatrix TriangleStressHelpers::generate_B() const {
    const double det_j = x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3);

    // From [1] 4.61
    BMatrix b;
    // clang-format off
    b << y(2, 3),     0.0, y(3, 1),     0.0, y(1, 2),     0.0,
             0.0, x(3, 2),     0.0, x(1, 3),     0.0, x(2, 1),
         x(3, 2), y(2, 3), x(1, 3), y(3, 1), x(2, 1), y(1, 2);
    // clang-format on
    return b / det_j;
}

//
// #############################################################################
//

double TriangleStressHelpers::area() const {
    // From [1] 4.70, I'm keeping the indices one-indexed like they do
    return 0.5f * abs(x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3));
}

//
// #############################################################################
//

GlobalKMatrix generate_global_stiffness_matrix(const SimulationContext &context) {
    // Generate the global stiffness matrix
    const size_t num_displacements = context.state.displacements.size();
    GlobalKMatrix K(num_displacements, num_displacements);

    // Each node can have up to 6 connections, this means we'll have at most a 6x6
    // matrix for each displacement
    K.reserve(6 * 6 * num_displacements);
    K.setZero();

    for (size_t i = 0; i < context.mesh.triangles.size(); ++i) {
        // No processing needed if the triangle is fixed
        if (context.cache.fixed_triangles[i]) continue;

        TriangleStressHelpers{context.mesh.triangles[i], context}.populate_local_stiffness(K);
    }
    return K;
}

//
// #############################################################################
//

SimulationContext::SimulationContext(common::Mesh mesh_) {
    mesh = std::move(mesh_);
    const size_t vertex_count = mesh.vertices.size();

    cache.vertex_to_displacements.resize(vertex_count, -1);
    size_t num_dynamic_displacements = 0;
    for (size_t i = 0; i < vertex_count; ++i) {
        if (mesh.fixed[i]) continue;

        cache.vertex_to_displacements[i] = num_dynamic_displacements;  // store before incrementing
        num_dynamic_displacements++;
    }

    state.displacements = Eigen::VectorXd::Zero(num_dynamic_displacements);
    state.velocities = Eigen::VectorXd::Zero(num_dynamic_displacements);
    state.accelerations = Eigen::VectorXd::Zero(num_dynamic_displacements);

    // Apply gravity (only to the Y coordinates)
    for (size_t i = 1; i < num_dynamic_displacements; i += 2) {
        state.accelerations[i] = -9.8;
    }

    cache.fixed_triangles.reserve(vertex_count);
    for (const auto &triangle : mesh.triangles) {
        // Fixed until proven otherwise
        bool fixed = true;

        for (size_t index : triangle.indices) {
            if (!mesh.fixed[index]) {
                // At least one index isn't fixed!
                fixed = false;
                break;
            }
        }

        cache.fixed_triangles.push_back(fixed);
    }

    cache.mass = generate_mass_matrix(num_dynamic_displacements, cache.vertex_to_displacements, mesh.mass);
    cache.damping = generate_damping_matrix(num_dynamic_displacements);

    GlobalKMatrix K = generate_global_stiffness_matrix(*this);

    // [2] Table 9.3 A.5
    K = K + MeshStepper::kA0 * cache.mass + MeshStepper::kA1 * cache.damping;

    // [2] Table 9.3 A.6
    cache.K_solver.compute(K);
    if (cache.K_solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition Failed!");
    }
}

//
// #############################################################################
//

State MeshStepper::step(SimulationContext &context) {
    const size_t num_displacements = context.state.displacements.size();

    // Copy terminology to match [2] Table 9.3
    const Eigen::MatrixXd &M = context.cache.mass;
    const Eigen::MatrixXd &C = context.cache.damping;
    const Eigen::VectorXd &U = context.state.displacements;
    const Eigen::VectorXd &U_vel = context.state.velocities;
    const Eigen::VectorXd &U_accel = context.state.accelerations;

    Eigen::VectorXd gravity = Eigen::VectorXd::Zero(num_displacements);
    for (size_t i = 1; i < num_displacements; i += 2) {
        gravity[i] = -9.8 * M(i, i);
    }

    // [2] Table 9.3 B.1
    Eigen::VectorXd R_hat =
        gravity + M * (kA0 * U + kA2 * U_vel + kA3 * U_accel) + C * (kA1 * U + kA4 * U_vel + kA5 * U_accel);

    // [2] Table 9.3 B.2
    Eigen::VectorXd U_next = context.cache.K_solver.solve(R_hat);
    if (context.cache.K_solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving Failed!");
    }

    // [2] Table 9.3 B.3
    Eigen::VectorXd U_accel_next = kA0 * (U_next - U) - kA2 * U_vel - kA3 * U_accel;
    Eigen::VectorXd U_vel_next = U_vel + kA6 * U_accel + kA7 * U_accel_next;

    // Now we can update our internal state!
    State next_state;
    next_state.displacements = std::move(U_next);
    next_state.velocities = std::move(U_vel_next);
    next_state.accelerations = std::move(U_accel_next);

    // To avoid numerical issues with things moving too quickly, we'll impose a
    // max velocity
    for (size_t i = 0; i < num_displacements; ++i) {
        double &velocity = next_state.velocities[i];

        constexpr double kTerminalVelocity = 100.0;
        velocity = std::clamp(velocity, -kTerminalVelocity, kTerminalVelocity);
    }

    return next_state;
}

//
// #############################################################################
//

std::vector<double> MeshStepper::evaluate_stresses() const { return {}; }

//
// #############################################################################
//

void Simulator::step(const double dt) {
    if (!context) return;

    const size_t num_steps = dt / MeshStepper::kDt;
    for (size_t i = 0; i < num_steps; ++i) {
        context->state = MeshStepper::step(*context);
    }

    destroy_stressful_triangles();
}

//
// #############################################################################
//

void Simulator::draw() const {
    if (!context) return;

    for (size_t i = 0; i < context->mesh.triangles.size(); ++i) {
        glBegin(GL_TRIANGLES);

        const common::Triangle &triangle = context->mesh.triangles[i];

        if (context->cache.fixed_triangles[i]) {
            glColor3f(0.f, 0.f, 1.f);
        } else {
            constexpr double kMaxStress = 5'000'000;
            const Eigen::Vector3d stress_v = compute_stress(triangle);
            const double stress = stress_v.norm();

            float red = static_cast<float>(stress / kMaxStress);  // 0 when stress is 0, 1 when stress is high
            float green =
                static_cast<float>((kMaxStress - stress) / kMaxStress);  // 1 when stress is 0, 0 when stress is high
            glColor3f(std::clamp(red, 0.f, 1.f), std::clamp(green, 0.f, 1.f), 0.0f);
        }

        // constexpr float kPxPerMeter = 200;
        constexpr float kPxPerMeter = static_cast<float>(common::kPxSize) / common::kBlockSize;
        const auto x = [&](uint8_t i) { return kPxPerMeter * context->get_coordinate(triangle.indices[i]); };
        const auto y = [&](uint8_t i) { return kPxPerMeter * context->get_coordinate(triangle.indices[i] + 1); };

        glVertex2f(x(0), y(0));
        glVertex2f(x(1), y(1));
        glVertex2f(x(2), y(2));
        glEnd();

        glBegin(GL_LINE_LOOP);
        glColor3f(0.1f, 0.1f, 0.1f);
        glVertex2f(x(0), y(0));
        glVertex2f(x(1), y(1));
        glVertex2f(x(2), y(2));
        glEnd();
    }
}

//
// #############################################################################
//

void Simulator::set_mesh(common::Mesh mesh) { context = std::make_unique<SimulationContext>(std::move(mesh)); }

//
// #############################################################################
//

void Simulator::destroy_stressful_triangles() {
    std::vector<double> stresses;
    stresses.reserve(context->mesh.triangles.size());
    double max_stress = 0;
    for (const auto &triangle : context->mesh.triangles) {
        double stress = compute_stress(triangle).norm();
        max_stress = std::max(max_stress, stress);
        stresses.emplace_back(stress);
    }

    constexpr double kMaxStress = 10'000'000;
    if (max_stress < kMaxStress) {
        return;
    }

    std::vector<common::Triangle> triangles;
    triangles.reserve(context->mesh.triangles.size());

    auto previous_state = context->state;
    for (size_t i = 0; i < context->mesh.triangles.size(); ++i) {
        if (stresses[i] < kMaxStress) {
            triangles.emplace_back(context->mesh.triangles[i]);
        }
    }
    context->mesh.triangles = std::move(triangles);

    // Using set_mesh will reset the displacements/velocities/accelerations so
    // we need to make sure to record them and reuse them for the next iteration
    set_mesh(std::move(context->mesh));
    context->state = std::move(previous_state);
}

//
// #############################################################################
//

Eigen::Vector3d Simulator::compute_stress(const common::Triangle &triangle) const {
    Eigen::Matrix<double, 6, 1> displacements;
    for (size_t i = 0; i < triangle.indices.size(); ++i) {
        displacements[2 * i] = context->get_displacement(triangle.indices[i]);
        displacements[2 * i + 1] = context->get_displacement(triangle.indices[i] + 1);
    }
    return generate_D() * TriangleStressHelpers{triangle, *context}.generate_B() * displacements;
}
