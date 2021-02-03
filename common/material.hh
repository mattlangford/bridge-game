#pragma once

#include <cstdint>
#include <string_view>
#include <tuple>

namespace common {

enum class Material : uint8_t {
    kNone = 0,
    kBrick = 1,
    kStone = 2,
    kRoad = 3,
};

struct Properties {
    /// Name of the material, like "road"
    std::string_view name;

    /// RGB color of the material
    std::tuple<float, float, float> color;

    /// Mass density of the material in units of kg/m^3
    double mass_density;

    /// Essentially the modulus of elasticity of the material, Measured in N/m^2, higher numbers mean stiffer materials
    double youngs_modulus;

    /// Describes how strains in one dimension cause strains in the other. Unitless
    double poissons_ratio;

    /// Maximum total stress before breaking. Measured in N/m^2.
    /// TODO: This should be split up between compressive and tensile stress
    double max_stress;

    /// Should this material be treated as fixed for the physics simulation
    bool fixed;
};

///
/// @brief Fetch the properties for the given material
///
const Properties& get_properties(const Material& m);
}  // namespace common
