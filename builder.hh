#include <optional>
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
        return 0.3;
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
        handler.add().key(GLFW_KEY_1, [this](GLFWwindow*, int) {
            drawing_cell_ = Cell::kBrick;
            data_[index(0, 0)] = drawing_cell_;
        });
        handler.add().key(GLFW_KEY_2, [this](GLFWwindow*, int) {
            drawing_cell_ = Cell::kRoad;
            data_[index(0, 0)] = drawing_cell_;
        });

        // Update cursor zoom
        handler.add().key(GLFW_KEY_Z, [this](GLFWwindow*, int) {
            cursor_zoom_++;
        });
        handler.add().key(GLFW_KEY_X, [this](GLFWwindow*, int) {
            cursor_zoom_ = 1;
        });

        handler.add().key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            erase_mode_ = true;
        });
        handler.add().release().key(GLFW_KEY_LEFT_CONTROL, [this](GLFWwindow*, int) {
            erase_mode_ = false;
        });

        handler.add().any_type().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_mouse_block(xpos, ypos);
        });

        // Draw
        handler.add().left_click([this](GLFWwindow*) {
            set_cell_at_mouse_block(erase_mode_ ? Cell::kNone : drawing_cell_);
        });
        handler.add().hold().any_modifier().move([this](GLFWwindow*, double xpos, double ypos) {
            set_cell_at_mouse_block(erase_mode_ ? Cell::kNone : drawing_cell_);
        });
    }

    Mesh convert_to_mesh() const
    {
        return {};
    }

private:
    size_t index(uint16_t w_block, uint16_t h_block) const
    {
        return w_block * num_h_blocks + h_block;
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
