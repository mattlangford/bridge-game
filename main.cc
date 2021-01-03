#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <GLFW/glfw3.h>

#include "builder.hh"
#include "events.hh"

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

    init_view();

    Builder builder;
    EventHandler handler;

    handler.add().key(GLFW_KEY_ESCAPE, [](GLFWwindow* window, int) { glfwSetWindowShouldClose(window, 1); });
    handler.add().key(GLFW_KEY_Z, [](GLFWwindow* window, int) { std::cout << "Zoom!\n"; });
    handler.add(EventState::kInit).move([](GLFWwindow*, double x, double) { std::cout << x << "x\n"; });
    handler.add(EventState::kInit).control().move([](GLFWwindow*, double, double y) { std::cout << y << "y\n"; });
    handler.add(EventState::kInit).right_click([](GLFWwindow*) { std::cout << "right click\n"; });
    handler.add(EventState::kInit).twice().right_click([](GLFWwindow*) { std::cout << "double right click\n"; });

    glfwSetMouseButtonCallback(window, route_mouse_button_callback);
    glfwSetCursorPosCallback(window, route_cursor_position_callback);
    glfwSetKeyCallback(window, route_key_callback);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        builder.draw();

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Terminate GLFW
    glfwTerminate();

    // Exit program
    exit(EXIT_SUCCESS);
}
