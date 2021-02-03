#include "common/mesh.hh"

namespace common {

//
// #############################################################################
//

std::vector<size_t> remove_orphaned_vertices(Mesh& mesh) {
    // Go through and figure out which vertices we can throw away
    std::vector<bool> to_keep(mesh.vertices.size(), false);
    size_t num_to_keep = 0;
    for (const auto& triangle : mesh.triangles) {
        for (size_t index : triangle.indices) {
            if (to_keep[index]) continue;

            to_keep[index] = true;      // x
            to_keep[index + 1] = true;  // y
            num_to_keep += 2;
        }
    }

    Mesh new_mesh;
    new_mesh.vertices.resize(num_to_keep, 0.0);
    new_mesh.fixed.resize(num_to_keep, false);
    new_mesh.mass.resize(num_to_keep, 0.0);

    // Generate the mapping between old vertices and new vertices and populate the vertex level mesh data
    std::vector<size_t> results(mesh.vertices.size(), -1);
    for (size_t from = 0, to = 0; from <= mesh.vertices.size(); ++from) {
        if (!to_keep[from]) {
            continue;
        }
        results[from] = to;

        new_mesh.vertices[to] = mesh.vertices[from];
        new_mesh.fixed[to] = mesh.fixed[from];
        new_mesh.mass[to] = mesh.mass[from];

        to++;
    }

    new_mesh.triangles = std::move(mesh.triangles);
    mesh = std::move(new_mesh);

    // Now adjust the triangle indices
    for (auto& triangle : mesh.triangles) {
        for (size_t& index : triangle.indices) {
            index = results[index];
        }
    }

    return results;
}

//
// #############################################################################
//

std::vector<size_t> remove_triangles(std::vector<size_t> triangle_indices, Mesh& mesh) {
    std::sort(triangle_indices.begin(), triangle_indices.end());
    for (auto i = triangle_indices.rbegin(); i != triangle_indices.rend(); ++i) {
        mesh.triangles.erase(mesh.triangles.begin() + *i);
    }

    return remove_orphaned_vertices(mesh);
}
}  // namespace common
