#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <optional>
#include <utility>
#include <vector>

#include <GLFW/glfw3.h>

static constexpr size_t kWidth = 1280;
static constexpr size_t kHeight = 720;

enum class Cell : uint8_t {
    kNone = 0,
    kBrick = 1
};

class Drawer {
public:
    static constexpr size_t kPxSize = 30;

    inline static Drawer* drawer;

public:
    Drawer()
    {
        data_.resize(num_w_blocks * num_h_blocks, Cell::kNone);
        drawer = this;
    }
    ~Drawer()
    {
        drawer = nullptr;
    }

    void draw() const
    {
        // Start by drawing all of the blocks that we have
        glColor3f(0.75f, 0.5f, 0.0f);
        glBegin(GL_TRIANGLES);
        for (uint16_t w_block = 0; w_block < num_w_blocks; w_block++) {
            for (uint16_t h_block = 0; h_block < num_h_blocks; h_block++) {
                if (data_[index(w_block, h_block)] == Cell::kNone) {
                    continue;
                }

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

        // Then draw the grid
        glBegin(GL_LINES);
        glColor3f(0.25f, 0.25f, 0.25f);
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

        // And if the mouse cursor is on the screen, we can draw that too
        if (mouse_block_) {
            if (erase_mode_)
                glColor4f(1.0f, 0.2f, 0.2f, 0.5f);
            else
                glColor4f(1.0f, 1.0f, 1.0f, 0.5f);

            const auto& [w_block, h_block] = *mouse_block_;
            uint16_t w_px = w_block * kPxSize;
            uint16_t h_px = h_block * kPxSize;

            uint16_t w_end = w_px + kPxSize;
            uint16_t h_end = h_px + kPxSize;

            glVertex2i(w_px, h_px); // top left
            glVertex2i(w_end, h_px); // top right

            glVertex2i(w_end, h_px); // top right
            glVertex2i(w_end, h_end); // bottom right

            glVertex2i(w_end, h_end); // bottom right
            glVertex2i(w_px, h_end); // bottom left

            glVertex2i(w_px, h_end); // bottom left
            glVertex2i(w_px, h_px); // top left
        }

        glEnd();
    }

    size_t index(uint16_t w_block, uint16_t h_block) const
    {
        return w_block * num_h_blocks + h_block;
    }

    void callback(GLFWwindow* window)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        if (xpos < 0 || ypos < 0 || xpos >= kWidth || ypos >= kHeight) {
            mouse_block_ = std::nullopt;
            return;
        }

        uint16_t w_px = static_cast<uint16_t>(xpos) / kPxSize;
        uint16_t h_px = static_cast<uint16_t>(kHeight - ypos) / kPxSize;
        mouse_block_ = std::make_pair(w_px, h_px);

        const bool control = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL);
        erase_mode_ = control;

        const bool pressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (!pressed) {
            return;
        }
        data_[index(w_px, h_px)] = erase_mode_ ? Cell::kNone : Cell::kBrick;
    }

private:
    const size_t num_w_blocks = kWidth / kPxSize;
    const size_t num_h_blocks = kHeight / kPxSize;

    bool erase_mode_ = false;
    std::optional<std::pair<uint16_t, uint16_t>> mouse_block_;

    // row major cell data
    std::vector<Cell> data_;
};

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE) {
        glfwSetWindowShouldClose(window, 1);
    }

    if (!Drawer::drawer) {
        return;
    }

    Drawer::drawer->callback(window);
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (!Drawer::drawer) {
        return;
    }

    Drawer::drawer->callback(window);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (!Drawer::drawer) {
        return;
    }

    Drawer::drawer->callback(window);
}

void init_view()
{
    // set up view
    // glViewport(0, 0, kWidth, kHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // this creates a canvas to do 2D drawing on
    glOrtho(0.0, kWidth, 0.0, kHeight, 0.0, 1.0);
}

int main(int argc, char* argv[])
{
    GLFWwindow* window;
    Drawer drawer;

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    window = glfwCreateWindow(kWidth, kHeight, "Window", NULL, NULL);
    if (!window) {
        std::cerr << "Unable to create window!\n";
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);

    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    init_view();

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.0f, 0.1f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        drawer.draw();

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Terminate GLFW
    glfwTerminate();

    // Exit program
    exit(EXIT_SUCCESS);
}
