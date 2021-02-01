#pragma once

#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "common/config.hh"
#include "common/material.hh"
#include "common/mesh.hh"
#include "renderer/events.hh"

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

void generate_ridge(BuildingContext &context);

BuildingContext generate_initial_building_context();

common::Mesh generate_mesh(const BuildingContext &context);

void draw_grid();

void draw(const DrawingContext &context);
void draw(const BuildingContext &context);

class Builder {
   public:
    Builder();

    void setup_callbacks(EventHandler &handler);
    void draw() const;

    const BuildingContext &building_context() const;

   private:
    void set_drawing_cell(common::Material cell);
    common::Material get_drawing_cell() const;

    void cursor_zoom_in();
    void cursor_zoom_reset();

    void set_erase_mode(bool enabled = true);
    bool get_erase_mode() const;

    void set_mouse_block(double xpos, double ypos);

    void set_cell_at_mouse_block(const common::Material &new_cell, bool force = false);

   private:
    DrawingContext drawing_context_;
    BuildingContext building_context_;
};
