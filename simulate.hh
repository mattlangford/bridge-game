#include <chrono>
#include <iostream>
#include <optional>
#include <vector>

#include "mesh.hh"

class Simlator {
public:
    void step()
    {
        if (vertex_velocities_.size() != mesh_.vertices.size()) {
            throw std::runtime_error("Invalid mesh!");
        }

        std::vector<bool> fixed_mask(mesh_.vertices.size(), false);
        for (Triangle& triangle : mesh_.triangles) {
            if (!triangle.metadata.fixed) {
                continue;
            }

            for (size_t index : triangle.indices) {
                fixed_mask[index] = true;
                fixed_mask[index + 1] = true;
            }
        }

        for (size_t i = 0; i < mesh_.vertices.size(); ++i) {
            if (i % 2 == 0)
                continue;
            if (fixed_mask[i])
                continue;

            vertex_velocities_[i] -= 0.1;
            mesh_.vertices[i] += vertex_velocities_[i];
        }
    }

    void draw() const
    {
        draw_mesh(mesh_);
    }

    void set_mesh(Mesh mesh)
    {
        mesh_ = std::move(mesh);
        vertex_velocities_.resize(mesh_.vertices.size());
        std::fill(vertex_velocities_.begin(), vertex_velocities_.end(), 0.0f);
    }

private:
    std::vector<float> vertex_velocities_;
    Mesh mesh_;
};
