#include "renderer/builder.hh"

#include <GLFW/glfw3.h>

#include "builder/context.hh"
#include "common/config.hh"

void draw_grid() {
    glColor3f(0.15f, 0.25f, 0.25f);
    glBegin(GL_LINES);
    for (size_t w_block = 0; w_block < common::kNumWBlocks; w_block++) {
        uint16_t w_px = w_block * common::kPxSize;
        glVertex2i(w_px, 0);
        glVertex2i(w_px, common::kHeight);
    }
    for (size_t h_block = 0; h_block < common::kNumHBlocks; h_block++) {
        uint16_t h_px = h_block * common::kPxSize;
        glVertex2i(0, h_px);
        glVertex2i(common::kWidth, h_px);
    }
    glEnd();
}

//
// #############################################################################
//

void draw(const DrawingContext &context) {
    if (!context.mouse_block) {
        return;
    }

    if (context.erase_mode)
        glColor3f(1.0f, 0.3f, 0.2f);
    else
        glColor3f(1.0f, 1.0f, 1.0f);

    glBegin(GL_LINES);
    const auto &[w_block, h_block] = *context.mouse_block;
    uint16_t w_px = w_block * common::kPxSize;
    uint16_t h_px = h_block * common::kPxSize;

    uint16_t w_end = w_px + context.cursor_zoom * common::kPxSize;
    uint16_t h_end = h_px + context.cursor_zoom * common::kPxSize;

    // Expand the edges just a bit
    w_px = w_px <= 0 ? w_px : w_px - 1;
    h_px = h_px <= 0 ? h_px : h_px - 1;
    w_end = w_end >= common::kWidth ? w_end : w_end + 1;
    h_end = h_end >= common::kHeight ? h_end : h_end + 1;

    glVertex2i(w_px, h_px);   // top left
    glVertex2i(w_end, h_px);  // top right

    glVertex2i(w_end, h_px);   // top right
    glVertex2i(w_end, h_end);  // bottom right

    glVertex2i(w_end, h_end);  // bottom right
    glVertex2i(w_px, h_end);   // bottom left

    glVertex2i(w_px, h_end);  // bottom left
    glVertex2i(w_px, h_px);   // top left
    glEnd();
}

//
// #############################################################################
//

void draw(const BuildingContext &context) {
    auto last_material = common::Material::kNone;

    glColor3f(0.75f, 0.5f, 0.0f);
    glBegin(GL_TRIANGLES);
    for (uint16_t w_block = 0; w_block < common::kNumWBlocks; w_block++) {
        for (uint16_t h_block = 0; h_block < common::kNumHBlocks; h_block++) {
            const common::Material material = context.data[context.index(w_block, h_block)];

            // No need to draw when it's none
            if (material == common::Material::kNone) {
                continue;
            }

            // Only update the drawing color if the material changes
            if (material != last_material) {
                const auto [r, g, b] = common::color(material);
                glColor3f(r, g, b);

                last_material = material;
            }

            uint16_t w_px = w_block * common::kPxSize;
            uint16_t h_px = h_block * common::kPxSize;

            uint16_t w_end = w_px + common::kPxSize;
            uint16_t h_end = h_px + common::kPxSize;

            // top half
            glVertex2i(w_px, h_px);   // top left
            glVertex2i(w_end, h_px);  // top right
            glVertex2i(w_px, h_end);  // bottom left

            // bottom half
            glVertex2i(w_end, h_end);  // bottom right
            glVertex2i(w_end, h_px);   // top right
            glVertex2i(w_px, h_end);   // bottom left
        }
    }
    glEnd();

    draw_grid();
}
