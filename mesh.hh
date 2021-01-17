#pragma once

#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>

//
// #############################################################################
//

struct Connection {
    size_t triangle_index;
    size_t point_index;
};

//
// #############################################################################
//

struct Triangle {
    static constexpr size_t kDim = 3;

    ///
    /// Which vertices belong to this triangle
    ///
    std::array<size_t, kDim> indices;

    ///
    /// External connections to other triangles for each vertex
    ///
    std::array<std::vector<Connection>, kDim> connections;
};

//
// #############################################################################
//

struct Mesh {
    ///
    /// Stores the vertices in a flat vector (x0, y0, x1, y1, ..., xn yn). The size will be 2 * 3 * triangles.size().
    ///
    std::vector<float> vertices;

    ///
    /// Groups vertices into triangles. Since the vertices vector is flat, this will only connect the first coordinates.
    ///
    std::vector<Triangle> triangles;
};

//
// #############################################################################
//

///
/// @brief Builds a mesh from integer coordinates
///
class MeshBuilder {

private:
    using Coordinate = uint32_t;
    using Coordinate2d = std::pair<Coordinate, Coordinate>;
    using Triangle = std::array<std::pair<Coordinate, Coordinate>, 3>;

public:
    MeshBuilder(size_t num_triangles_guess = 0)
    {
        triangles_.reserve(num_triangles_guess);
    }

public:
    void add_triangle(const Coordinate2d& c0, const Coordinate2d& c1, const Coordinate2d& c2)
    {
        triangles_.emplace_back(Triangle { c0, c1, c2 });
    }

    Mesh finalize() const
    {
        // Init the final mesh, reserving all the space we'll need
        Mesh mesh;
        mesh.vertices.reserve(3 * triangles_.size());
        mesh.triangles.reserve(triangles_.size());

        populate_triangles(mesh);
        populate_connections(mesh);

        return mesh;
    }

private:
    void populate_triangles(Mesh& mesh) const
    {
        // We'll use this to establish connections between vertices, this works since we're using integer coordinates
        std::unordered_map<Coordinate2d, size_t, Coordinate2dHash> visited;

        // We'll create the triangles in order they were added. It's possible we can do better here
        for (const auto& triangle : triangles_) {
            auto& mesh_triangle = mesh.triangles.emplace_back();
            for (size_t i = 0; i < triangle.size(); ++i) {
                const Coordinate2d& coord = triangle[i];
                const size_t next_vertex_index = mesh.vertices.size();

                // Attempt to emplace this coordinate in the visited map
                auto [it, emplaced] = visited.emplace(std::make_pair(coord, next_vertex_index));
                if (emplaced) {
                    // If this coordinates hadn't been seen before we can add it to the vertices vector
                    mesh.vertices.emplace_back(coord.first);
                    mesh.vertices.emplace_back(coord.second);
                }

                mesh_triangle.indices[i] = it->second;
            }
        }
    }

    void populate_connections(Mesh& mesh) const
    {
        // Each vertex may be connected to different triangles. In order to figure out what is connected to what, we'll
        // go through a first pass and generate a bunch of these ConnectionHelpers then go through and populate the
        // actual data in the mesh.
        struct ConnectionHelper {
            // All of the connections associated with this entry
            std::vector<Connection> connections;

            // Each of the above connections will go into each of these here
            std::vector<std::vector<Connection>*> destinations;
        };

        // Since we only store connections to the first coordinate in each pair, we only need half the size here
        // TODO(mlangford): This can probably be made better and avoid a crap ton of copies
        std::vector<ConnectionHelper> connections(mesh.vertices.size() / 2);
        for (size_t triangle_index = 0; triangle_index < mesh.triangles.size(); ++triangle_index) {
            auto& triangle = mesh.triangles[triangle_index];
            for (size_t coord_index = 0; coord_index < triangle.indices.size(); ++coord_index) {
                const size_t vertex_index = triangle.indices[coord_index];
                auto& connection = triangle.connections[coord_index];

                auto& helper = connections[vertex_index / 2];
                helper.connections.emplace_back(Connection { triangle_index, coord_index });
                helper.destinations.emplace_back(&connection);
            }
        }

        // Now push all of the connections into all of the destinations in each helper
        for (const ConnectionHelper& helper : connections) {
            for (size_t destination_index = 0; destination_index < helper.destinations.size(); ++destination_index) {
                auto& destination = *helper.destinations[destination_index];
                for (size_t connection_index = 0; connection_index < helper.connections.size(); ++connection_index) {
                    // Don't put the same connection into it's own destination
                    if (connection_index == destination_index) {
                        continue;
                    }

                    const auto& connection = helper.connections[connection_index];
                    destination.push_back(connection);
                }
            }
        }
    }

private:
    struct Coordinate2dHash {
        size_t operator()(Coordinate2d input) const
        {
            size_t result = 0;
            result |= input.first;
            result |= static_cast<size_t>(input.second) << sizeof(Coordinate);
            return result;
        }
    };

    ///
    /// We'll store them as triangles as we're building, but then convert to a flat vector at the end
    ///
    std::vector<Triangle> triangles_;
};

//
// #############################################################################
//
