#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <GLFW/glfw3.h>

#include "builder.hh"
#include "events.hh"
#include "mesh.hh"

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

    glfwSetMouseButtonCallback(window, route_mouse_button_callback);
    glfwSetCursorPosCallback(window, route_cursor_position_callback);
    glfwSetKeyCallback(window, route_key_callback);

    MeshBuilder mesh_builder;
    mesh_builder.add_triangle({ 300, 300 }, { 400, 300 }, { 300, 400 }, { false, 1.0 });
    mesh_builder.add_triangle({ 400, 300 }, { 300, 400 }, { 400, 400 }, { false, 0.5 });
    mesh_builder.add_triangle({ 600, 600 }, { 632, 550 }, { 700, 100 }, { false, 0.1 });
    const auto mesh = mesh_builder.finalize();

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // builder.draw();
        draw_mesh(mesh);

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Terminate GLFW
    glfwTerminate();

    // Exit program
    exit(EXIT_SUCCESS);
}
