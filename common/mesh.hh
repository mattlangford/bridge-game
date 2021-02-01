#pragma once

#include <array>
#include <compare>
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
} // namespace common
