#pragma once
#include <GLFW/glfw3.h>

struct DrawingContext;
struct BuildingContext;

void draw_grid();
void draw(const DrawingContext &context);
void draw(const BuildingContext &context);
