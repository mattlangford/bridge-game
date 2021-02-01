#include "common/material.hh"

namespace common {

std::tuple<float, float, float> color(const Material &c) {
    switch (c) {
        case Material::kBrick:
            return {0.75f, 0.5f, 0.0f};
        case Material::kStone:
            return {0.4f, 0.4f, 0.5f};
        case Material::kRoad:
            return {0.1f, 0.1f, 0.1f};
        case Material::kNone:
        default:
            return {0.5f, 0.5f, 0.5f};
    }
}

//
// #############################################################################
//

double mass(const Material &c) {
    switch (c) {
        case Material::kBrick:
            // Assume 50 bricks per m^3 and 3.1kg per brick
            return 50.0 * 3.1;
        // TODO: The rest of these
        case Material::kStone:
            return 0.9;
        case Material::kRoad:
            return 0.7;
        case Material::kNone:
        default:
            return 0.0;
    }
}

//
// #############################################################################
//

double youngs_modulus(const Material &) { return 3.7 * 1E5; }

//
// #############################################################################
//

double poissons_ration(const Material &) { return 0.1; }

//
// #############################################################################
//

bool fixed(const Material &c) { return c == Material::kStone; }
}  // namespace common
