#pragma once

#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "builder/mesh_builder.hh"
#include "common/mesh.hh"
#include "renderer/events.hh"

static constexpr size_t kWidth = 1280;
static constexpr size_t kHeight = 720;

enum class Cell : uint8_t {
    kNone = 0,
    kBrick = 1,
    kStone = 2,
    kRoad = 3,
};

void set_color_from_cell(const Cell& c)
{
    switch (c) {
    case Cell::kBrick:
        glColor3f(0.75f, 0.5f, 0.0f);
        return;
    case Cell::kStone:
        glColor3f(0.4f, 0.4f, 0.5f);
        return;
    case Cell::kRoad:
        glColor3f(0.1f, 0.1f, 0.1f);
        return;
    case Cell::kNone:
    default:
        return;
    }
}

double get_mass_from_cell(const Cell& c)
{
    switch (c) {
    case Cell::kBrick:
        // Assume 50 bricks per m^3 and 3.1kg per brick
        return 50.0 * 3.1;
    case Cell::kStone:
        return 0.9;
    case Cell::kRoad:
        return 0.7;
    case Cell::kNone:
    default:
        return 0.0;
    }
}

struct DrawingContext {
    /// Allows toggling between erase mode and draw mode
    bool erase_mode = false;

    /// How zoomed in the cursor is
    size_t cursor_zoom = 1;

    /// If the mouse is on screen, where is it at
    std::optional<std::pair<uint16_t, uint16_t>> mouse_block;

    /// What is the current cell being drawn
    Cell drawing_cell = Cell::kBrick;
};

struct BuildingContext {
    /// How big each block is for rendering
    static constexpr size_t kPxSize = 15; // px
    /// How big each block is for simulation
    static constexpr size_t kBlockSize = 5; // meters

    /// How many blocks there are on the screen
    static constexpr size_t kNumWBlocks = kWidth / kPxSize;
    static constexpr size_t kNumHBlocks = kHeight / kPxSize;

    /// Row major cell data starting from the bottom left
    std::vector<Cell> data;

    /// Helper functions to make indexing easier
    size_t index(uint16_t w_block, uint16_t h_block) const
    {
        return w_block * kNumHBlocks + h_block;
    }
    std::pair<uint16_t, uint16_t> reverse_index(size_t index) const
    {
        auto [q, r] = std::ldiv(index, kNumHBlocks);
        return { q, r };
    }
};

void generate_ridge(BuildingContext& context)
{
    // The bottom of the ridge
    for (size_t w_block = 0; w_block < BuildingContext::kNumWBlocks; w_block++) {
        auto index = [&](size_t row) { return context.index(w_block, row); };
        context.data[index(0)] = Cell::kStone;
        context.data[index(1)] = Cell::kStone;
        context.data[index(BuildingContext::kNumHBlocks / 2)] = Cell::kBrick;
    }

    // On the sides
    for (size_t h_block = 0; h_block < BuildingContext::kNumHBlocks / 2; h_block++) {
        auto index = [&](size_t col) { return context.index(col, h_block); };
        context.data[index(0)] = Cell::kStone;
        context.data[index(1)] = Cell::kStone;
        context.data[index(BuildingContext::kNumWBlocks - 1)] = Cell::kStone;
        context.data[index(BuildingContext::kNumWBlocks - 2)] = Cell::kStone;
    }
}

BuildingContext generate_initial_building_context()
{
    BuildingContext context;
    context.data.resize(BuildingContext::kNumHBlocks * BuildingContext::kNumHBlocks, Cell::kNone);
    generate_ridge(context);
    return context;
}

Mesh generate_mesh(const BuildingContext& context)
{
    MeshBuilder builder;

    std::unordered_set<size_t> visited;

    // Find a place to start a new mesh, we'll consider everything that's touching to be in the same mesh
    for (size_t index = 0; index < context.data.size(); ++index) {
        // Simple breadth first search
        std::queue<size_t> bfs;

        auto add_index = [&](const size_t index) {
            if (index >= context.data.size())
                return false; // Don't worry if it's off screen
            if (context.data[index] == Cell::kNone)
                return false; // Don't worry if it's not populated
            auto [it, emplaced] = visited.emplace(index);
            if (emplaced)
                bfs.push(index); // If it wasn't in the visited list we can queue it up
            return emplaced;
        };

        // If this first one fails we've either already visited this node or it's none
        if (!add_index(index)) {
            continue;
        }

        std::vector<size_t> indices;
        while (!bfs.empty()) {
            const size_t this_index = bfs.front();
            bfs.pop();

            indices.push_back(this_index);

            // Right neighbor
            add_index(this_index + 1);
            // Upper neighbor
            add_index(this_index + BuildingContext::kNumHBlocks);
            // Left neighbor
            add_index(this_index - 1);
            // Bottom neighbor
            add_index(this_index - BuildingContext::kNumHBlocks);
        }

        // Now we can go through and add the triangles to the meshes
        for (const size_t this_index : indices) {
            const Cell& cell = context.data[this_index];

            Metadata metadata;
            metadata.fixed = cell == Cell::kStone;

            constexpr auto kBlockSize = BuildingContext::kBlockSize;
            constexpr float kVolume = kBlockSize * kBlockSize * kBlockSize / 2.0;

            metadata.mass = kVolume * get_mass_from_cell(cell);

            const auto [w_block, h_block] = context.reverse_index(this_index);

            const uint16_t w_m = w_block * kBlockSize;
            const uint16_t h_m = h_block * kBlockSize;

            const uint16_t w_end = w_m + kBlockSize;
            const uint16_t h_end = h_m + kBlockSize;

            // top half
            builder.add_triangle({ w_m, h_m }, { w_end, h_m }, { w_m, h_end }, metadata);

            // bottom half
            builder.add_triangle({ w_end, h_end }, { w_end, h_m }, { w_m, h_end }, metadata);
        }
    }

    return builder.finalize();
}

void draw_grid()
{
    glColor3f(0.15f, 0.25f, 0.25f);
    glBegin(GL_LINES);
    for (size_t w_block = 0; w_block < BuildingContext::kNumWBlocks; w_block++) {
        uint16_t w_px = w_block * BuildingContext::kPxSize;
        glVertex2i(w_px, 0);
        glVertex2i(w_px, kHeight);
    }
    for (size_t h_block = 0; h_block < BuildingContext::kNumHBlocks; h_block++) {
        uint16_t h_px = h_block * BuildingContext::kPxSize;
        glVertex2i(0, h_px);
        glVertex2i(kWidth, h_px);
    }
    glEnd();
}

void draw(const DrawingContext& context)
{
    if (!context.mouse_block) {
        return;
    }

    if (context.erase_mode)
        glColor3f(1.0f, 0.3f, 0.2f);
    else
        glColor3f(1.0f, 1.0f, 1.0f);

    glBegin(GL_LINES);
    const auto& [w_block, h_block] = *context.mouse_block;
    uint16_t w_px = w_block * BuildingContext::kPxSize;
    uint16_t h_px = h_block * BuildingContext::kPxSize;

    uint16_t w_end = w_px + context.cursor_zoom * BuildingContext::kPxSize;
    uint16_t h_end = h_px + context.cursor_zoom * BuildingContext::kPxSize;

    // Expand the edges just a bit
    w_px = w_px <= 0 ? w_px : w_px - 1;
    h_px = h_px <= 0 ? h_px : h_px - 1;
    w_end = w_end >= kWidth ? w_end : w_end + 1;
    h_end = h_end >= kHeight ? h_end : h_end + 1;

    glVertex2i(w_px, h_px); // top left
    glVertex2i(w_end, h_px); // top right

    glVertex2i(w_end, h_px); // top right
    glVertex2i(w_end, h_end); // bottom right

    glVertex2i(w_end, h_end); // bottom right
    glVertex2i(w_px, h_end); // bottom left

    glVertex2i(w_px, h_end); // bottom left
    glVertex2i(w_px, h_px); // top left
    glEnd();
}

void draw(const BuildingContext& context)
{
    glColor3f(0.75f, 0.5f, 0.0f);
    glBegin(GL_TRIANGLES);
    for (uint16_t w_block = 0; w_block < BuildingContext::kNumWBlocks; w_block++) {
        for (uint16_t h_block = 0; h_block < BuildingContext::kNumHBlocks; h_block++) {
            const Cell cell = context.data[context.index(w_block, h_block)];
            if (cell == Cell::kNone) {
                continue;
            }
            set_color_from_cell(cell);

            uint16_t w_px = w_block * BuildingContext::kPxSize;
            uint16_t h_px = h_block * BuildingContext::kPxSize;

            uint16_t w_end = w_px + BuildingContext::kPxSize;
            uint16_t h_end = h_px + BuildingContext::kPxSize;

            // top half
            glVertex2i(w_px, h_px); // top left
            glVertex2i(w_end, h_px); // top right
            glVertex2i(w_px, h_end); // bottom left

            // bottom half
            glVertex2i(w_end, h_end); // bottom right
            glVertex2i(w_end, h_px); // top right
            glVertex2i(w_px, h_end); // bottom left
        }
    }
    glEnd();

    draw_grid();
}

class Builder {
public:
    Builder()
        : building_context_(generate_initial_building_context())
    {
    }

    const BuildingContext& building_context() const { return building_context_; }
    void draw() const
    {
        ::draw(building_context_);
        ::draw(drawing_context_);
    }

    void setup_callbacks(EventHandler& handler)
    {
        // Updating the drawing mode
        handler.add(EventState::kBuild).key(GLFW_KEY_1, [this](GLFWwindow*, int) {
            set_drawing_cell(Cell::kBrick);
        });
        handler.add(EventState::kBuild).key(GLFW_KEY_2, [this](GLFWwindow*, int) {
            set_drawing_cell(Cell::kRoad);
        });

        // Update cursor zoom
        handler.add(EventState::kBuild).key(GLFW_KEY_Z, [this](GLFWwindow*, int) {
            cursor_zoom_in();
        });
        handler.add(EventState::kBuild).key(GLFW_KEY_X, [this](GLFWwindow*, int) {
            cursor_zoom_reset();
        });

        handler.add(EventState::kBuild).key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            set_erase_mode();
        });
        handler.add(EventState::kBuild).release().key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            set_erase_mode(false);
        });

        handler.add(EventState::kBuild).any_type().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_mouse_block(xpos, ypos);
        });

        // Draw
        handler.add(EventState::kBuild).any_modifier().left_click([this](GLFWwindow*) {
            set_cell_at_mouse_block(get_drawing_cell());
        });
        handler.add(EventState::kBuild).hold().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_cell_at_mouse_block(get_drawing_cell());
        });
    }

private:
    void set_drawing_cell(Cell cell)
    {
        drawing_context_.drawing_cell = cell;
    }

    Cell get_drawing_cell() const
    {
        return get_erase_mode() ? Cell::kNone : drawing_context_.drawing_cell;
    }

    void cursor_zoom_in()
    {
        drawing_context_.cursor_zoom++;
    }
    void cursor_zoom_reset()
    {
        drawing_context_.cursor_zoom = 1;
    }

    void set_erase_mode(bool enabled = true)
    {
        drawing_context_.erase_mode = true;
    }

    bool get_erase_mode() const
    {
        return drawing_context_.erase_mode;
    }

    void set_mouse_block(double xpos, double ypos)
    {
        // Don't worry about anything else here if the mouse is off screen
        if (xpos < 0 || ypos < 0 || xpos >= kWidth || ypos >= kHeight) {
            drawing_context_.mouse_block = std::nullopt;
            return;
        }

        uint16_t w_block = static_cast<uint16_t>(xpos) / BuildingContext::kPxSize;
        uint16_t h_block = static_cast<uint16_t>(kHeight - ypos) / BuildingContext::kPxSize;
        drawing_context_.mouse_block = std::make_pair(w_block, h_block);
    }

    void set_cell_at_mouse_block(const Cell& new_cell, bool force = false)
    {
        if (!drawing_context_.mouse_block)
            return;

        const auto& [w_block, h_block] = *drawing_context_.mouse_block;
        for (size_t w = 0; w < drawing_context_.cursor_zoom; ++w) {
            for (size_t h = 0; h < drawing_context_.cursor_zoom; ++h) {
                Cell& cell = building_context_.data[building_context_.index(w_block + w, h_block + h)];
                if (cell != Cell::kStone || force) {
                    cell = new_cell;
                }
            }
        }
    }

private:
    DrawingContext drawing_context_;
    BuildingContext building_context_;
};
