#include <optional>
#include <utility>
#include <vector>

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

    // Lazy single callback for all actions (regardless of if they're keyboard or mouse)
    void callback(GLFWwindow* window)
    {
        // Update if we're erasing or not
        const bool control = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL);
        erase_mode_ = control;

        // Updating the drawing mode
        if (glfwGetKey(window, GLFW_KEY_1)) {
            drawing_cell_ = Cell::kBrick;
        } else if (glfwGetKey(window, GLFW_KEY_2)) {
            drawing_cell_ = Cell::kRoad;
        }
        data_[index(0, 0)] = drawing_cell_;

        // Update cursor zoom
        if (glfwGetKey(window, GLFW_KEY_Z)) {
            cursor_zoom_++;
        } else if (glfwGetKey(window, GLFW_KEY_X)) {
            cursor_zoom_ = 1;
        }

        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        // Don't worry about anything else here if the mouse is off screen
        if (xpos < 0 || ypos < 0 || xpos >= kWidth || ypos >= kHeight) {
            mouse_block_ = std::nullopt;
            return;
        }

        uint16_t w_block = static_cast<uint16_t>(xpos) / kPxSize;
        uint16_t h_block = static_cast<uint16_t>(kHeight - ypos) / kPxSize;
        mouse_block_ = std::make_pair(w_block, h_block);

        // Make sure we're pressed before we start drawing
        const bool pressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (!pressed) {
            return;
        }

        // Populate a cell. Stone can't be overwritten
        const auto set_cell = [this](uint16_t w_block, uint16_t h_block) {
            Cell& cell = data_[index(w_block, h_block)];
            if (cell != Cell::kStone) {
                cell = erase_mode_ ? Cell::kNone : drawing_cell_;
            }
        };

        // Since the cursor may be zoomed, deal with that
        for (size_t w = 0; w < cursor_zoom_; ++w) {
            for (size_t h = 0; h < cursor_zoom_; ++h) {
                set_cell(w_block + w, h_block + h);
            }
        }
    }

private:
    size_t index(uint16_t w_block, uint16_t h_block) const
    {
        return w_block * num_h_blocks + h_block;
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
