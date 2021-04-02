#pragma once

#include "common/material.hh"

namespace common {
/// Canvas size
static constexpr size_t kWidth = 1280;
static constexpr size_t kHeight = 720;

/// Each simulation update will progress forward by this amount
static constexpr double kSimulationDt = 1.0 / 240.0;
/// Every frame will be this far apart in time (by calling the simulation update multiple times)
static constexpr double kRenderDt = 1.0 / 60.0;

/// How big each block is for rendering
static constexpr size_t kPxSize = 20;  // px
/// How big each block is for simulation
static constexpr size_t kBlockSize = 2;  // meters

/// How many blocks there are on the screen
static constexpr size_t kNumWBlocks = kWidth / kPxSize;
static constexpr size_t kNumHBlocks = kHeight / kPxSize;

/// Triangle destruction is a bit buggy at the moment, so only enable it for testing
static constexpr bool kEnableTriangleDestruction = false;

/// Used for the damping matrix
static constexpr double kDampingFactor = 0.0;

/// Max speed for falling objects, used to fix numerical issues with very rapidly moving objects
static constexpr double kTerminalVelocity = 100.0;

/// Definitions for material properties
inline Properties get_brick_properties() {
    Properties prop;
    prop.name = "brick";
    prop.color = {0.75f, 0.5f, 0.0f};
    // Assume 50 bricks per m^3 and 3.1kg per brick
    prop.mass_density = 50.0 * 3.1;
    prop.youngs_modulus = 3.7 * 1E7;
    prop.poissons_ratio = 0.1;
    prop.max_stress = 6'00'000;
    prop.fixed = false;
    return prop;
}
inline Properties get_stone_properties() {
    Properties prop;
    prop.name = "stone";
    prop.color = {0.4f, 0.4f, 0.5f};
    prop.fixed = true;

    // Generally these should be unused
    prop.mass_density = 0.0;
    prop.youngs_modulus = 0.0;
    prop.poissons_ratio = 0.0;
    prop.max_stress = 1E10;
    return prop;
}
inline Properties get_road_properties() {
    Properties prop;
    prop.name = "road";
    prop.color = {0.1f, 0.1f, 0.1f};
    // Copying the values for prop for now
    prop.mass_density = 50.0 * 3.1;
    prop.youngs_modulus = 3.7 * 1E8;
    prop.poissons_ratio = 0.1;
    prop.max_stress = 5'000'000;
    prop.fixed = false;
    return prop;
}
}  // namespace common
