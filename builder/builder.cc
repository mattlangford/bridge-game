#include "builder/builder.hh"

#include "builder/mesh_builder.hh"

void generate_ridge(BuildingContext &context) {
    // The bottom of the ridge
    for (size_t w_block = 0; w_block < common::kNumWBlocks; w_block++) {
        auto index = [&](size_t row) { return context.index(w_block, row); };
        context.data[index(0)] = common::Material::kStone;
        context.data[index(common::kNumHBlocks / 2)] = common::Material::kBrick;
    }

    // On the sides
    for (size_t h_block = 0; h_block < common::kNumHBlocks / 2; h_block++) {
        auto index = [&](size_t col) { return context.index(col, h_block); };
        context.data[index(0)] = common::Material::kStone;
        context.data[index(common::kNumWBlocks - 1)] = common::Material::kStone;
    }
}

//
// #############################################################################
//

BuildingContext generate_initial_building_context() {
    BuildingContext context;
    context.data.resize(common::kNumHBlocks * common::kNumWBlocks, common::Material::kNone);
    generate_ridge(context);
    return context;
}

//
// #############################################################################
//

common::Mesh generate_mesh(const BuildingContext &context) {
    MeshBuilder builder;

    std::unordered_set<size_t> visited;

    // Find a place to start a new mesh, we'll consider everything that's touching
    // to be in the same mesh
    for (size_t index = 0; index < context.data.size(); ++index) {
        // Simple breadth first search
        std::queue<size_t> bfs;

        auto add_index = [&](const size_t index) {
            if (index >= context.data.size()) return false;                    // Don't worry if it's off screen
            if (context.data[index] == common::Material::kNone) return false;  // Don't worry if it's not populated
            auto [it, emplaced] = visited.emplace(index);
            if (emplaced) bfs.push(index);  // If it wasn't in the visited list we can queue it up
            return emplaced;
        };

        // If this first one fails we've either visited this node or it's kNone
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
            add_index(this_index + common::kNumHBlocks);
            // Left neighbor
            add_index(this_index - 1);
            // Bottom neighbor
            add_index(this_index - common::kNumHBlocks);
        }

        // Now we can go through and add the triangles to the meshes
        for (const size_t this_index : indices) {
            const common::Material &material = context.data[this_index];

            const auto [w_block, h_block] = context.reverse_index(this_index);

            const uint16_t w_m = w_block * common::kBlockSize;
            const uint16_t h_m = h_block * common::kBlockSize;

            const uint16_t w_end = w_m + common::kBlockSize;
            const uint16_t h_end = h_m + common::kBlockSize;

            // top half
            builder.add_triangle({w_m, h_m}, {w_end, h_m}, {w_m, h_end}, material);

            // bottom half
            builder.add_triangle({w_end, h_end}, {w_end, h_m}, {w_m, h_end}, material);
        }
    }

    return builder.finalize();
}

//
// #############################################################################
//

Builder::Builder() : building_context_(generate_initial_building_context()) {}

//
// #############################################################################
//

const BuildingContext &Builder::building_context() const { return building_context_; }
const DrawingContext &Builder::drawing_context() const { return drawing_context_; }

//
// #############################################################################
//

void Builder::setup_callbacks(EventHandler &handler) {
    // Updating the drawing mode
    handler.add(EventState::kBuild).key(GLFW_KEY_1, [this](GLFWwindow *, int) {
        set_drawing_cell(common::Material::kBrick);
    });
    handler.add(EventState::kBuild).key(GLFW_KEY_2, [this](GLFWwindow *, int) {
        set_drawing_cell(common::Material::kRoad);
    });

    // Update cursor zoom
    handler.add(EventState::kBuild).key(GLFW_KEY_Z, [this](GLFWwindow *, int) { cursor_zoom_in(); });
    handler.add(EventState::kBuild).key(GLFW_KEY_X, [this](GLFWwindow *, int) { cursor_zoom_reset(); });

    handler.add(EventState::kBuild).key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow *, int) { set_erase_mode(); });
    handler.add(EventState::kBuild).release().key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow *, int) {
        set_erase_mode(false);
    });

    handler.add(EventState::kBuild).any_type().any_modifier().move([this](GLFWwindow *, double xpos, double ypos) {
        set_mouse_block(xpos, ypos);
    });

    // Draw
    handler.add(EventState::kBuild).any_modifier().left_click([this](GLFWwindow *) {
        set_cell_at_mouse_block(get_drawing_cell());
    });
    handler.add(EventState::kBuild).hold().any_modifier().move([this](GLFWwindow *, double, double) {
        set_cell_at_mouse_block(get_drawing_cell());
    });
}

void Builder::set_drawing_cell(common::Material cell) { drawing_context_.drawing_cell = cell; }

common::Material Builder::get_drawing_cell() const {
    return get_erase_mode() ? common::Material::kNone : drawing_context_.drawing_cell;
}

void Builder::cursor_zoom_in() { drawing_context_.cursor_zoom++; }
void Builder::cursor_zoom_reset() { drawing_context_.cursor_zoom = 1; }

void Builder::set_erase_mode(bool enabled) { drawing_context_.erase_mode = enabled; }

bool Builder::get_erase_mode() const { return drawing_context_.erase_mode; }

void Builder::set_mouse_block(double xpos, double ypos) {
    // Don't worry about anything else here if the mouse is off screen
    if (xpos < 0 || ypos < 0 || xpos >= common::kWidth || ypos >= common::kHeight) {
        drawing_context_.mouse_block = std::nullopt;
        return;
    }

    uint16_t w_block = static_cast<uint16_t>(xpos) / common::kPxSize;
    uint16_t h_block = static_cast<uint16_t>(common::kHeight - ypos) / common::kPxSize;
    drawing_context_.mouse_block = std::make_pair(w_block, h_block);
}

void Builder::set_cell_at_mouse_block(const common::Material &new_cell, bool force) {
    if (!drawing_context_.mouse_block) return;

    const auto &[w_block, h_block] = *drawing_context_.mouse_block;
    for (size_t w = 0; w < drawing_context_.cursor_zoom; ++w) {
        for (size_t h = 0; h < drawing_context_.cursor_zoom; ++h) {
            if (common::kPxSize * (w_block + w) >= common::kWidth) continue;
            if (common::kPxSize * (h_block + h) >= common::kHeight) continue;
            common::Material &cell = building_context_.data[building_context_.index(w_block + w, h_block + h)];
            if (cell != common::Material::kStone || force) {
                cell = new_cell;
            }
        }
    }
}
