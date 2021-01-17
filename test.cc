#include "mesh.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#define FORMAT(lhs, rhs, name)                                                          \
    "Line: " << __LINE__ << ", " << name << " failed since these are not equal:\n" \
                                      << "  " #lhs << " = " << lhs << "\n"              \
                                      << "  " #rhs << " = " << rhs

#define ASSERT_EQ(lhs, rhs)                  \
    if ((lhs) != (rhs)) {                    \
        std::stringstream ss;                \
        ss << FORMAT(lhs, rhs, "ASSERT_EQ"); \
        throw std::runtime_error(ss.str());  \
    }
#define EXPECT_EQ(lhs, rhs)                                      \
    if ((lhs) != (rhs)) {                                        \
        std::cout << FORMAT(lhs, rhs, "EXPECT_EQ") << std::endl; \
    }

void test_mesh()
{
    MeshBuilder builder;
    builder.add_triangle({ 0, 0 }, { 1, 0 }, { 0, 1 });
    builder.add_triangle({ 1, 0 }, { 0, 1 }, { 1, 1 });
    builder.add_triangle({ 1, 0 }, { 1, 1 }, { 2, 1 });

    const auto mesh = builder.finalize();

    // Vertex check
    {
        std::vector<size_t> expected { 0, 0, 1, 0, 0, 1, 1, 1, 2, 1 };
        ASSERT_EQ(mesh.vertices.size(), expected.size());
        for (size_t i = 0; i < mesh.vertices.size(); ++i)
            EXPECT_EQ(mesh.vertices[i], expected[i]);
    }

    // Triangle vertex index check
    {
        std::array<size_t, 3> tri1 { 0, 2, 4 };
        std::array<size_t, 3> tri2 { 2, 4, 6 };
        std::array<size_t, 3> tri3 { 2, 6, 8 };
        ASSERT_EQ(mesh.triangles.size(), 3);
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_EQ(mesh.triangles[0].indices[i], tri1[i]);
            EXPECT_EQ(mesh.triangles[1].indices[i], tri2[i]);
            EXPECT_EQ(mesh.triangles[2].indices[i], tri3[i]);
        }
    }

    // Triangle connection index check
    {
        ASSERT_EQ(mesh.triangles.size(), 3);

        auto& tri0 = mesh.triangles[0];
        ASSERT_EQ(tri0.connections[0].size(), 0);
        ASSERT_EQ(tri0.connections[1].size(), 2);
        ASSERT_EQ(tri0.connections[2].size(), 1);

        auto& tri1 = mesh.triangles[1];
        ASSERT_EQ(tri1.connections[0].size(), 2);
        ASSERT_EQ(tri1.connections[1].size(), 1);
        ASSERT_EQ(tri1.connections[2].size(), 1);

        auto& tri2 = mesh.triangles[2];
        ASSERT_EQ(tri2.connections[0].size(), 2);
        ASSERT_EQ(tri2.connections[1].size(), 1);
        ASSERT_EQ(tri2.connections[2].size(), 0);

        // spot check of a few of them
        ASSERT_EQ(tri0.connections[2][0].triangle_index, 1);
        ASSERT_EQ(tri0.connections[2][0].point_index, 1);
        ASSERT_EQ(tri1.connections[0][0].triangle_index, 0);
        ASSERT_EQ(tri1.connections[0][0].point_index, 1);
        ASSERT_EQ(tri2.connections[1][0].triangle_index, 1);
        ASSERT_EQ(tri2.connections[1][0].point_index, 2);
    }
}

int main(int argc, char* argv[])
{
    std::unordered_map<std::string, std::function<void(void)>> tests;
    tests["mesh"] = &test_mesh;

    if (argc == 2) {
        const std::string name = argv[1];
        auto it = tests.find(name);
        if (it == std::end(tests)) {
            std::cout << "No tests found named '" << name << "'!\n";
            return 1;
        }

        std::cout << "Test: '" << it->first << "'\n";
        try {
            it->second();
        } catch (const std::exception& ex) {
            std::cout << "FAILED: '" << ex.what() << "'\n";
            return 1;
        }
        std::cout << "PASSED\n";
    } else {
        std::cout << "Running all tests!\n";
        for (const auto& [name, func] : tests) {
            std::cout << "Test: '" << name << "'\n";
            try {
                func();
            } catch (const std::exception& ex) {
                std::cout << "FAILED: '" << ex.what() << "'\n";
                continue;
            }
            std::cout << "PASSED\n";
        }
    }

    return 0;
}
