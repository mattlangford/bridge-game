#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "events.hh"
#include "mesh.hh"

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

class Builder {
public:
    static constexpr size_t kPxSize = 10;

    inline static Builder* drawer;

public:
    Builder()
    {
        data_.resize(num_w_blocks * num_h_blocks, Cell::kNone);
        generate_ridge();

        drawer = this;
    }
    ~Builder()
    {
        drawer = nullptr;
    }

    void draw() const
    {
        draw_cells();
        draw_grid();
        draw_mouse();
    }

    void setup_callbacks(EventHandler& handler)
    {
        // Updating the drawing mode
        handler.add(EventState::kBuild).key(GLFW_KEY_1, [this](GLFWwindow*, int) {
            drawing_cell_ = Cell::kBrick;
            data_[index(0, 0)] = drawing_cell_;
        });
        handler.add(EventState::kBuild).key(GLFW_KEY_2, [this](GLFWwindow*, int) {
            drawing_cell_ = Cell::kRoad;
            data_[index(0, 0)] = drawing_cell_;
        });
        // data_[index(0, 0)] = drawing_cell_;

        // Update cursor zoom
        handler.add(EventState::kBuild).key(GLFW_KEY_Z, [this](GLFWwindow*, int) {
            cursor_zoom_++;
        });
        handler.add(EventState::kBuild).key(GLFW_KEY_X, [this](GLFWwindow*, int) {
            cursor_zoom_ = 1;
        });

        handler.add(EventState::kBuild).key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            erase_mode_ = true;
        });
        handler.add(EventState::kBuild).release().key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            erase_mode_ = false;
        });

        handler.add(EventState::kBuild).any_type().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_mouse_block(xpos, ypos);
        });

        // Draw
        handler.add(EventState::kBuild).left_click([this](GLFWwindow*) {
            set_cell_at_mouse_block(erase_mode_ ? Cell::kNone : drawing_cell_);
        });
        handler.add(EventState::kBuild).hold().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_cell_at_mouse_block(erase_mode_ ? Cell::kNone : drawing_cell_);
        });
    }

    Mesh generate_mesh() const
    {
        MeshBuilder builder;

        std::unordered_set<size_t> visited;

        // Find a place to start a new mesh, we'll consider everything that's touching to be in the same mesh
        for (size_t index = 0; index < data_.size(); ++index) {
            // Simple breadth first search
            std::queue<size_t> bfs;

            auto add_index = [&](const size_t index) {
                if (index >= data_.size())
                    return false; // Don't worry if it's off screen
                if (data_[index] == Cell::kNone)
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
                add_index(this_index + num_h_blocks);
                // Left neighbor
                add_index(this_index - 1);
                // Bottom neighbor
                add_index(this_index - num_h_blocks);
            }

            // Now we can go through and add the triangles to the meshes
            for (const size_t this_index : indices) {
                const Cell& cell = data_[this_index];

                constexpr size_t kBlockSize = 1; // meters

                Metadata metadata;
                metadata.fixed = cell == Cell::kStone;

                constexpr float kVolume = kBlockSize * kBlockSize / 2.0;
                metadata.mass = kVolume * get_mass_from_cell(cell);

                const auto [w_block, h_block] = reverse_index(this_index);

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

private:
    size_t index(uint16_t w_block, uint16_t h_block) const
    {
        return w_block * num_h_blocks + h_block;
    }
    std::pair<uint16_t, uint16_t> reverse_index(size_t index) const
    {
        auto [q, r] = std::ldiv(index, num_h_blocks);
        return { q, r };
    }

    void set_mouse_block(double xpos, double ypos)
    {
        // Don't worry about anything else here if the mouse is off screen
        if (xpos < 0 || ypos < 0 || xpos >= kWidth || ypos >= kHeight) {
            mouse_block_ = std::nullopt;
            return;
        }

        uint16_t w_block = static_cast<uint16_t>(xpos) / kPxSize;
        uint16_t h_block = static_cast<uint16_t>(kHeight - ypos) / kPxSize;
        mouse_block_ = std::make_pair(w_block, h_block);
    }

    void set_cell_at_mouse_block(const Cell& new_cell, bool force = false)
    {
        if (!mouse_block_)
            return;

        const auto& [w_block, h_block] = *mouse_block_;
        for (size_t w = 0; w < cursor_zoom_; ++w) {
            for (size_t h = 0; h < cursor_zoom_; ++h) {
                Cell& cell = data_[index(w_block + w, h_block + h)];
                if (cell != Cell::kStone || force) {
                    cell = new_cell;
                }
            }
        }
    }

    void draw_cells() const
    {
        glColor3f(0.75f, 0.5f, 0.0f);
        glBegin(GL_TRIANGLES);
        for (uint16_t w_block = 0; w_block < num_w_blocks; w_block++) {
            for (uint16_t h_block = 0; h_block < num_h_blocks; h_block++) {
                const Cell cell = data_[index(w_block, h_block)];
                if (cell == Cell::kNone) {
                    continue;
                }
                set_color_from_cell(cell);

                uint16_t w_px = w_block * kPxSize;
                uint16_t h_px = h_block * kPxSize;

                uint16_t w_end = w_px + kPxSize;
                uint16_t h_end = h_px + kPxSize;

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
    }

    void draw_grid() const
    {
        glColor3f(0.15f, 0.25f, 0.25f);
        glBegin(GL_LINES);
        for (size_t w_block = 0; w_block < num_w_blocks; w_block++) {
            uint16_t w_px = w_block * kPxSize;
            glVertex2i(w_px, 0);
            glVertex2i(w_px, kHeight);
        }
        for (size_t h_block = 0; h_block < num_h_blocks; h_block++) {
            uint16_t h_px = h_block * kPxSize;
            glVertex2i(0, h_px);
            glVertex2i(kWidth, h_px);
        }
        glEnd();
    }

    void draw_mouse() const
    {
        if (!mouse_block_) {
            return;
        }

        if (erase_mode_)
            glColor3f(1.0f, 0.3f, 0.2f);
        else
            glColor3f(1.0f, 1.0f, 1.0f);

        glBegin(GL_LINES);
        const auto& [w_block, h_block] = *mouse_block_;
        uint16_t w_px = w_block * kPxSize;
        uint16_t h_px = h_block * kPxSize;

        uint16_t w_end = w_px + cursor_zoom_ * kPxSize;
        uint16_t h_end = h_px + cursor_zoom_ * kPxSize;

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

    void generate_ridge()
    {
        data_[500] = Cell::kStone;
        return;
        // The bottom of the ridge
        for (size_t w_block = 0; w_block < num_w_blocks; w_block++) {
            data_[index(w_block, 0)] = Cell::kStone;
            data_[index(w_block, 1)] = Cell::kStone;
        }

        for (size_t h_block = 0; h_block < num_h_blocks / 2; h_block++) {
            data_[index(0, h_block)] = Cell::kStone;
            data_[index(1, h_block)] = Cell::kStone;
            data_[index(num_w_blocks - 1, h_block)] = Cell::kStone;
            data_[index(num_w_blocks - 2, h_block)] = Cell::kStone;
        }
    }

private:
    const size_t num_w_blocks = kWidth / kPxSize;
    const size_t num_h_blocks = kHeight / kPxSize;

    bool erase_mode_ = false;
    size_t cursor_zoom_ = 1;
    std::optional<std::pair<uint16_t, uint16_t>> mouse_block_;
    Cell drawing_cell_ = Cell::kBrick;

    // Row major cell data starting from the bottom left
    std::vector<Cell> data_;
};
