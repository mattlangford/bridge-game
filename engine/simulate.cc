#include "engine/simulate.hh"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

#include "common/config.hh"
#include "iterate.hh"

static Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

//
// #############################################################################
//

Eigen::MatrixXd generate_mass_matrix(size_t num_displacements, const std::vector<size_t> &vertex_to_displacements,
                                     const std::vector<double> &mass) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(num_displacements, num_displacements);
    M.setZero();
    for (auto [index, mass] : it::zip(vertex_to_displacements, mass)) {
        if (index < num_displacements) {
            M(index, index) += mass;
        }
    }
    return M;
}

//
// #############################################################################
//

Eigen::MatrixXd generate_damping_matrix(size_t vertex_count) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(vertex_count, vertex_count);
    C.diagonal().fill(common::kDampingFactor);
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
    // Don't include the displacements, effectively calculating the initial area
    auto x = [&](size_t i, size_t j){
        auto i_coord = context.mesh.vertices[triangle.indices[i - 1]];
        auto j_coord = context.mesh.vertices[triangle.indices[j - 1]];
        return i_coord - j_coord;
    };
    auto y = [&](size_t i, size_t j){
        auto i_coord = context.mesh.vertices[triangle.indices[i] + 1];
        auto j_coord = context.mesh.vertices[triangle.indices[j] + 1];
        return i_coord - j_coord;
    };
    // From [1] 4.70, I'm keeping the indices one-indexed like they do
    return 0.5f * abs(x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3));
}

//
// #############################################################################
//

DMatrix TriangleStressHelpers::generate_D() const {
    const auto &properties = common::get_properties(triangle.material);
    const double E = properties.youngs_modulus;
    const double v = properties.poissons_ratio;

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
    for (auto [index, fixed] : it::zip(cache.vertex_to_displacements, mesh.fixed)) {
        if (fixed) continue;

        index = num_dynamic_displacements;  // store before incrementing
        num_dynamic_displacements++;
    }

    state.num_nodes = num_dynamic_displacements;
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

    cache.K = generate_global_stiffness_matrix(*this);

    // [2] Table 9.3 A.5/6
    cache.K_solver.compute(cache.K + MeshStepper::kA0 * cache.mass + MeshStepper::kA1 * cache.damping);
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
    next_state.num_nodes = static_cast<size_t>(next_state.displacements.size());

    // To avoid numerical issues with things moving too quickly, we'll impose a
    // max velocity
    for (size_t i = 0; i < next_state.num_nodes; ++i) {
        double &velocity = next_state.velocities[i];
        velocity = std::clamp(velocity, -common::kTerminalVelocity, common::kTerminalVelocity);
    }

    return next_state;
}

//
// #############################################################################
//

void Simulator::step(const double dt) {
    if (!context) return;

    const size_t num_steps = dt / MeshStepper::kDt;
    for (size_t i = 0; i < num_steps; ++i) {
        context->state = MeshStepper::step(*context);
    }

    if (common::kEnableTriangleDestruction) {
        destroy_stressful_triangles();
    }

    // Populate the stresses at the very end
    const auto &triangles = context->mesh.triangles;
    auto &triangle_stresses = context->cache.triangle_stresses;
    triangle_stresses.resize(triangles.size(), 3);
    triangle_stresses.setZero();
    for (auto [i, triangle] : it::enumerate(triangles)) {
        triangle_stresses.row(i) = compute_stress(triangle);
    }

    // Compute total energy, this should always be close to 0
    // const auto &state = context->state;
    // const auto &mass = context->cache.mass;
    // const double kinetic = 0.5 * state.velocities.transpose() * mass * state.velocities;
    // const double strain = 0.5 * state.displacements.transpose() * context->cache.K * state.displacements;
    // double gravity = 0;
    // for (size_t i = 1; i < state.num_nodes; i += 2) {
    //     gravity += -9.8 * mass.diagonal()[i] * -state.displacements[i];
    // }
    // const double total_energy = kinetic + strain + gravity;
    // std::cout << "Total energy : " << total_energy << " (KE:" << kinetic << ", gPE: " << gravity << ", sPE: " <<
    // strain << ")\n";
}

//
// #############################################################################
//

const SimulationContext *Simulator::simulation_context() const { return context.get(); }

//
// #############################################################################
//

void Simulator::set_mesh(common::Mesh mesh) { context = std::make_unique<SimulationContext>(std::move(mesh)); }

//
// #############################################################################
//

void Simulator::destroy_stressful_triangles() {
    std::vector<size_t> too_stressful;
    for (auto [i, triangle] : it::enumerate(context->mesh.triangles)) {
        double stress = compute_stress(triangle).norm();
        if (stress > common::get_properties(triangle.material).max_stress) {
            too_stressful.push_back(i);
        }
    }

    // TODO: This would be nice, but then it gets weird how we call set_mesh()
    // if (too_stressful.empty()) {
    //     set_mesh(std::move(context->mesh));
    //     return;
    // }

    // Using set_mesh will reset the displacements/velocities/accelerations so
    // we need to make sure to record them and reuse them for the next iteration
    auto previous_state = std::move(context->state);
    auto previous_vertex_to_displacements = std::move(context->cache.vertex_to_displacements);

    auto mapping = common::remove_triangles(too_stressful, context->mesh);
    set_mesh(std::move(context->mesh));

    for (auto [from, to] : it::enumerate(mapping)) {
        // Was this vertex kept? If not we don't have to worry about it
        if (to >= context->mesh.vertices.size()) {
            continue;
        }

        // Was this a dynamic vertex from the previous cycle? If not there will be no displacements
        const size_t dynamic_from = previous_vertex_to_displacements[from];
        if (dynamic_from >= previous_state.num_nodes) {
            continue;
        }

        // Is the kept vertex a dynamic index? If not there is nowhere to put it. If this passed the "dynamic_from"
        // check then it should generally pass this one as well.
        size_t dynamic_to = context->cache.vertex_to_displacements[to];
        if (dynamic_to >= context->state.num_nodes) {
            continue;
        }

        context->state.displacements[dynamic_to] = previous_state.displacements[dynamic_from];
        context->state.velocities[dynamic_to] = previous_state.velocities[dynamic_from];
        context->state.accelerations[dynamic_to] = previous_state.accelerations[dynamic_from];
    }
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
    TriangleStressHelpers helper{triangle, *context};
    return helper.generate_D() * helper.generate_B() * displacements;
}
