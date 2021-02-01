#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "common/config.hh"
#include "common/material.hh"

struct DrawingContext {
    /// Allows toggling between erase mode and draw mode
    bool erase_mode = false;

    /// How zoomed in the cursor is
    size_t cursor_zoom = 1;

    /// If the mouse is on screen, where is it at
    std::optional<std::pair<uint16_t, uint16_t>> mouse_block;

    /// What is the current cell being drawn
    common::Material drawing_cell = common::Material::kBrick;
};

struct BuildingContext {
    /// Row major cell data starting from the bottom left
    std::vector<common::Material> data;

    /// Helper functions to make indexing easier
    inline size_t index(uint16_t w_block, uint16_t h_block) const { return w_block * common::kNumHBlocks + h_block; }
    inline std::pair<uint16_t, uint16_t> reverse_index(size_t index) const {
        auto [q, r] = std::ldiv(index, common::kNumHBlocks);
        return {q, r};
    }
};
