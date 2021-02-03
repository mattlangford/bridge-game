#include "common/material.hh"

#include <stdexcept>

#include "common/config.hh"

namespace common {

static Properties kBrick = get_brick_properties();
static Properties kStone = get_stone_properties();
static Properties kRoad = get_road_properties();

//
// #############################################################################
//

const Properties& get_properties(const Material& m) {
    switch (m) {
        case (Material::kBrick):
            return kBrick;
        case (Material::kStone):
            return kStone;
        case (Material::kRoad):
            return kRoad;
        default:
            throw std::runtime_error("Invalid material! Cannot get properties.");
    }
}
}  // namespace common
