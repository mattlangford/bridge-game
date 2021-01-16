#include "mesh.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#define ASSERT_EQ(lhs, rhs)                                                                          \
    if ((lhs) != (rhs)) {                                                                            \
        std::stringstream ss;                                                                        \
        ss << "Line: " << __LINE__ << ", ASSERT_EQ failed since lhs: " << lhs << " != rhs: " << rhs; \
        throw std::runtime_error(ss.str());                                                          \
    }
#define EXPECT_EQ(lhs, rhs)                                                                                         \
    if ((lhs) != (rhs)) {                                                                                           \
        std::cout << "Line: " << __LINE__ << ", ASSERT_EQ failed since lhs: " << lhs << " != rhs: " << rhs << "\n"; \
    }

void test_mesh()
{
    MeshBuilder builder;
    builder.add_triangle({ 0, 0 }, { 1, 0 }, { 0, 1 });
    builder.add_triangle({ 1, 0 }, { 0, 1 }, { 1, 1 });

    const auto mesh = builder.finalize();

    // Vertex check
    {
        std::vector<size_t> expected { 0, 0, 1, 0, 0, 1, 1, 1 };
        ASSERT_EQ(mesh.vertices.size(), 2 * 4);
        for (size_t i = 0; i < mesh.vertices.size(); ++i)
            ASSERT_EQ(mesh.vertices[i], expected[i]);
    }

    // Triangle vertex index check
    {
        std::array<size_t, 3> tri1 { 0, 2, 4 };
        std::array<size_t, 3> tri2 { 2, 4, 6 };
        ASSERT_EQ(mesh.triangles.size(), 2);
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_EQ(mesh.triangles.front().indices[i], tri1[i]);
            EXPECT_EQ(mesh.triangles.back().indices[i], tri2[i]);
        }
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
