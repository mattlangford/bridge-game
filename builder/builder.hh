#pragma once

#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "builder/context.hh"
#include "common/config.hh"
#include "common/material.hh"
#include "common/mesh.hh"
#include "renderer/events.hh"

void generate_ridge(BuildingContext &context);

BuildingContext generate_initial_building_context();

common::Mesh generate_mesh(const BuildingContext &context);

void draw_grid();

class Builder {
   public:
    Builder();

    void setup_callbacks(EventHandler &handler);

    const BuildingContext &building_context() const;
    const DrawingContext &drawing_context() const;

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
