#pragma once

#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "common/material.hh"

namespace common {
//
// #############################################################################
//

struct Triangle {
    ///
    /// Which vertices belong to this triangle
    ///
    static constexpr size_t kDim = 3;
    std::array<size_t, kDim> indices;

    Material material;
};

//
// #############################################################################
//

struct Mesh {
    ///
    /// Stores the vertices in a flat vector (x0, y0, x1, y1, ..., xn yn). The
    /// size will be 2 * 3 * triangles.size().
    ///
    std::vector<float> vertices;

    ///
    /// Groups vertices into triangles. Since the vertices vector is flat, this
    /// will only connect the first coordinates.
    ///
    std::vector<Triangle> triangles;

    ///
    /// Stores if the vertex at the same index is fixed or not. Fixed vertices
    /// always have zeroed offsets in physics
    ///
    std::vector<bool> fixed;

    ///
    /// Stores each vertex's mass (in kg/m3)
    ///
    std::vector<double> mass;
};

///
/// @brief Remove vertices that have no triangles associated with them.
///
/// @returns A mapping from old to new indices, where the ith element is the new vertex index (or -1 if removed)
///
std::vector<size_t> remove_orphaned_vertices(Mesh& mesh);

///
/// @brief Remove triangles if their index is provided in the triangle indices vector. This also removes vertices.
/// NOTE: This may reorder the triangles
///
/// @returns A mapping from old to new indices, where the ith element is the new vertex index (or -1 if removed)
///
std::vector<size_t> remove_triangles(std::vector<size_t> triangle_indices, Mesh& mesh);
}  // namespace common
