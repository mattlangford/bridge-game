#pragma once

#include <cstdint>
#include <tuple>

namespace common {

enum class Material : uint8_t {
    kNone = 0,
    kBrick = 1,
    kStone = 2,
    kRoad = 3,
};

///
/// @brief Get (r,g,b) for the color of the material. Each entry will be between 0.0f and 1.0f
///
std::tuple<float, float, float> color(const Material& c);

///
/// @brief Get mass in units of kg/m^3
///
double mass(const Material& c);

///
/// @brief Get youngs modulus for the material in N/M^2
///
double youngs_modulus(const Material& c);

///
/// @brief Get poissons ratio (unitless)
///
double poissons_ration(const Material& c);

///
/// @brief Should we consider this material to be immoveable
///
bool fixed(const Material& c);
}
