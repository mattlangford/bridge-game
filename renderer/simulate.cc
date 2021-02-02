#include "renderer/simulate.hh"

#include <GLFW/glfw3.h>

#include "common/config.hh"
#include "engine/context.hh"

void draw(const SimulationContext& context) {
    const auto& cache = context.cache;
    const auto& mesh = context.mesh;

    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const common::Triangle& triangle = mesh.triangles[i];

        glBegin(GL_TRIANGLES);

        if (cache.fixed_triangles[i]) {
            glColor3f(0.f, 0.f, 1.f);
        } else {
            const double stress = cache.triangle_stresses.row(i).norm();

            constexpr float kMaxStress = 10'000;
            float red = static_cast<float>(stress / kMaxStress);  // 0 when stress is 0, 1 when stress is high
            float green =
                static_cast<float>((kMaxStress - stress) / kMaxStress);  // 1 when stress is 0, 0 when stress is high
            glColor3f(std::clamp(red, 0.f, 1.f), std::clamp(green, 0.f, 1.f), 0.0f);
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
