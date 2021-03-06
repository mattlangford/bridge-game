#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "builder/builder.hh"
#include "common/config.hh"
#include "engine/simulate.hh"
#include "renderer/builder.hh"
#include "renderer/events.hh"
#include "renderer/simulate.hh"

void init_view() {
    // set up view
    // glViewport(0, 0, kWidth, kHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // this creates a canvas to do 2D drawing on
    glOrtho(0.0, common::kWidth, 0.0, common::kHeight, 0.0, 1.0);
}

void print_state(EventState state) {
    std::cout << "EventHandler::set_state() to '" << event_state_to_string(state) << "'\n";
}

using Clock = std::chrono::high_resolution_clock;
void print_fps() {
    // Meh, global variables here just to make it easy
    static size_t frame_counter = 0;
    static Clock::time_point last_time = Clock::now();

    static constexpr size_t kFPSFrames = 100;
    if (++frame_counter % kFPSFrames == 0) {
        auto now = Clock::now();

        const double fps = kFPSFrames / std::chrono::duration<double>(now - last_time).count();
        std::cout << "Average fps: " << fps << "\n";

        last_time = now;
    }
}

int main() {
    EventHandler handler;

    GLFWwindow *window;

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    window = glfwCreateWindow(common::kWidth, common::kHeight, "Window", NULL, NULL);
    if (!window) {
        std::cerr << "Unable to create window!\n";
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    init_view();

    Builder builder;
    Simulator simulator;

    // Set up all the callbacks, note that some of these capture by reference but
    // since the EventHandler outlives everything in the main function, this
    // should be fine.
    handler.add_state_callback(&print_state);
    handler.add().key(GLFW_KEY_ESCAPE, [](GLFWwindow *window, int) { glfwSetWindowShouldClose(window, 1); });
    handler.add().key(GLFW_KEY_SPACE, [&](GLFWwindow *, int) {
        auto state = handler.get_state();
        switch (state) {
            case EventState::kInit:
            case EventState::kSimulate:
                handler.set_state(EventState::kBuild);
                break;
            case EventState::kBuild:
                handler.set_state(EventState::kSimulate);
                break;
        }
    });
    handler.add_state_callback(EventState::kSimulate,
                               [&](EventState) { simulator.set_mesh(generate_mesh(builder.building_context())); });
    builder.setup_callbacks(handler);

    bool step = true;
    handler.add().key(GLFW_KEY_P, [&](GLFWwindow *, int) { step = !step; });

    glfwSetMouseButtonCallback(window, route_mouse_button_callback);
    glfwSetCursorPosCallback(window, route_cursor_position_callback);
    glfwSetKeyCallback(window, route_key_callback);

    handler.set_state(EventState::kBuild);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        auto state = handler.get_state();
        switch (state) {
            case EventState::kBuild: {
                draw(builder.building_context());
                draw(builder.drawing_context());
                break;
            }
            case EventState::kSimulate: {
                if (step) simulator.step(common::kRenderDt);
                // step = false;
                if (auto ptr = simulator.simulation_context()) draw(*ptr);
                break;
            }
            default:
                continue;
        }

        glfwSwapBuffers(window);
        print_fps();
        glfwPollEvents();
    }

    // Terminate GLFW
    glfwTerminate();

    // Exit program
    exit(EXIT_SUCCESS);
}
