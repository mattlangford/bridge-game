#include "renderer/simulate.hh"

#include <GLFW/glfw3.h>

#include "common/config.hh"
#include "common/material.hh"
#include "engine/context.hh"
#include "engine/simulate.hh"

void draw_k_matrix(const SimulationContext& context)
{
    const auto& cache = context.cache;
    const auto& mesh = context.mesh;

    double max = context.cache.K.coeffs().maxCoeff();
    glBegin(GL_LINES);
    std::cout << Eigen::MatrixXd(context.cache.K) << "\n";
    for (size_t from = 0; from < mesh.vertices.size(); from+=2)
    {
        auto d_from = cache.vertex_to_displacements[from];
        if (d_from >= context.state.num_nodes) continue;
        for (size_t to = 0; to < mesh.vertices.size(); to+=2)
        {
            auto d_to = cache.vertex_to_displacements[to];
            if (d_to >= context.state.num_nodes) continue;

            if (d_from == d_to) continue;
            const double element = abs(context.cache.K.coeff(d_from, d_to));
            if (element < 1E-3) continue;

            constexpr float kPxPerMeter = static_cast<float>(common::kPxSize) / common::kBlockSize;
            float start_x = kPxPerMeter * context.get_coordinate(from);
            float start_y = kPxPerMeter * context.get_coordinate(from + 1);
            float end_x = kPxPerMeter * context.get_coordinate(to);
            float end_y = kPxPerMeter * context.get_coordinate(to + 1);

            float color = element / max;
            glColor3f(color, color, color);
            std::cout << "from: " << d_from << " to: " << d_to << " start: " << start_x << ", " << start_y << " to end: " << end_x << ", " <<  end_y << " mag: " << element <<"\n";
            glVertex2f(start_x, start_y);
            glVertex2f(end_x, end_y);
        }
    }
    glEnd();
    std::cout << "==============\n";
}

void draw(const SimulationContext& context) {
    const auto& cache = context.cache;
    const auto& mesh = context.mesh;

    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const common::Triangle& triangle = mesh.triangles[i];

        glBegin(GL_TRIANGLES);

        if (cache.fixed_triangles[i]) {
            glColor3f(0.f, 0.f, 1.f);
        } else {
            TriangleStressHelpers helper{triangle, context};
            const double zero_stress = 0.5 * common::kBlockSize * common::kBlockSize;
            const double stress = abs(helper.area() - zero_stress); //cache.triangle_stresses.row(i).norm();

            // Here we scale the actual max stress down a bit so it turns bright red before breaking
            const float max_stress = 0.02; // common::get_properties(triangle.material).max_stress;
            if (stress > max_stress)
            {
                glColor3f(0.f, 0.f, 0.f);
            }
            else
            {
                // 0 when stress is 0, 1 when stress is high
                const float red = static_cast<float>(stress / max_stress);
                // 1 when stress is 0, 0 when stress is high
                const float green = static_cast<float>((max_stress - stress) / max_stress);
                glColor3f(std::clamp(red, 0.f, 1.f), std::clamp(green, 0.f, 1.f), 0.0f);
            }
        }

        constexpr float kPxPerMeter = static_cast<float>(common::kPxSize) / common::kBlockSize;
        const auto x = [&](uint8_t i) { return kPxPerMeter * context.get_coordinate(triangle.indices[i]); };
        const auto y = [&](uint8_t i) { return kPxPerMeter * context.get_coordinate(triangle.indices[i] + 1); };

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
