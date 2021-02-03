#pragma once

#include <vector>

#include "builder/builder.hh"
#include "common/mesh.hh"

///
/// @brief Builds a mesh from integer coordinates
///
class MeshBuilder {
   private:
    using Coordinate = uint32_t;
    using Coordinate2d = std::pair<Coordinate, Coordinate>;
    struct BuildingTriangle {
        std::array<std::pair<Coordinate, Coordinate>, 3> coords;
        common::Material material;
    };

   public:
    MeshBuilder(size_t num_triangles_guess = 0) { triangles_.reserve(num_triangles_guess); }

   public:
    void add_triangle(Coordinate2d c0, Coordinate2d c1, Coordinate2d c2,
                      common::Material material = common::Material::kNone) {
        triangles_.emplace_back(BuildingTriangle{{c0, c1, c2}, material});
    }

    common::Mesh finalize() const {
        // Init the final mesh, reserving all the space we'll need
        common::Mesh mesh;
        mesh.vertices.reserve(3 * triangles_.size());
        mesh.triangles.reserve(triangles_.size());
        mesh.fixed.reserve(triangles_.size());
        mesh.mass.reserve(triangles_.size());

        populate_triangles(mesh);
        populate_fixed(mesh);
        populate_mass(mesh);

        return mesh;
    }

   private:
    void populate_triangles(common::Mesh &mesh) const {
        // We'll use this to establish connections between vertices, this works
        // since we're using integer coordinates
        std::unordered_map<Coordinate2d, size_t, Coordinate2dHash> visited;

        // We'll create the triangles in order they were added. It's possible we can
        // do better here
        for (const BuildingTriangle &triangle : triangles_) {
            auto &mesh_triangle = mesh.triangles.emplace_back();
            mesh_triangle.material = triangle.material;

            for (size_t i = 0; i < triangle.coords.size(); ++i) {
                const Coordinate2d &coord = triangle.coords[i];
                const size_t next_vertex_index = mesh.vertices.size();

                // Attempt to emplace this coordinate in the visited map
                auto [it, emplaced] = visited.emplace(std::make_pair(coord, next_vertex_index));
                if (emplaced) {
                    // If this coordinates hadn't been seen before we can add it to the
                    // vertices vector
                    mesh.vertices.emplace_back(coord.first);
                    mesh.vertices.emplace_back(coord.second);
                }

                mesh_triangle.indices[i] = it->second;
            }
        }
    }

    void populate_fixed(common::Mesh &mesh) const {
        mesh.fixed.resize(mesh.vertices.size(), false);
        for (size_t i = 0; i < triangles_.size(); ++i) {
            const auto &material = triangles_[i].material;
            if (material != common::Material::kStone) {
                continue;
            }

            for (size_t index : mesh.triangles[i].indices) {
                mesh.fixed[index] = true;      // x
                mesh.fixed[index + 1] = true;  // y
            }
        }
    }

    void populate_mass(common::Mesh &mesh) const {
        mesh.mass.resize(mesh.vertices.size(), 0.0);

        constexpr auto kBlockSize = common::kBlockSize;
        constexpr float kVolume = kBlockSize * kBlockSize * kBlockSize / 2.0;

        // Naively split the mass between each vertex. NOTE there are 6 vertices per
        // triangle (x1, y1, x2, y2, x3, y3) This could probably be improved by
        // integrating over the area of the triangles.
        for (size_t i = 0; i < triangles_.size(); ++i) {
            const auto sixth_mass = kVolume * common::get_properties(triangles_[i].material).mass_density / 6.0;
            for (size_t index : mesh.triangles[i].indices) {
                mesh.mass[index] += sixth_mass;
                mesh.mass[index + 1] += sixth_mass;
            }
        }
    }

   private:
    struct Coordinate2dHash {
        size_t operator()(Coordinate2d input) const {
            size_t result = 0;
            result |= input.first;
            result |= static_cast<size_t>(input.second) << sizeof(Coordinate);
            return result;
        }
    };

    ///
    /// We'll store them as triangles as we're building, but then convert to a
    /// flat vector at the end
    ///
    std::vector<BuildingTriangle> triangles_;
};
