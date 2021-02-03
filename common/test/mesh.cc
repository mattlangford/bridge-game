#include "common/mesh.hh"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "builder/mesh_builder.hh"

namespace common {

TEST(Mesh, remove_orphaned_vertices) {
    MeshBuilder builder;
    builder.add_triangle({0, 0}, {1, 1}, {2, 2}, Material::kBrick);
    builder.add_triangle({1, 1}, {2, 2}, {3, 3}, Material::kStone);

    auto mesh = builder.finalize();

    EXPECT_EQ(mesh.vertices.size(), 8);
    EXPECT_EQ(mesh.fixed.size(), 8);
    EXPECT_EQ(mesh.mass.size(), 8);

    // Make sure nothing was removed
    auto mapping = remove_orphaned_vertices(mesh);
    EXPECT_EQ(mesh.vertices.size(), 8);
    EXPECT_EQ(mesh.fixed.size(), 8);
    EXPECT_EQ(mesh.mass.size(), 8);
    EXPECT_THAT(mapping, testing::ElementsAre(0, 1, 2, 3, 4, 5, 6, 7));

    // Remove the first triangle
    ASSERT_EQ(mesh.triangles.size(), 2);
    mesh.triangles.erase(mesh.triangles.begin());
    ASSERT_EQ(mesh.triangles.size(), 1);

    mapping = remove_orphaned_vertices(mesh);
    EXPECT_EQ(mesh.vertices.size(), 6);
    EXPECT_EQ(mesh.fixed.size(), 6);
    EXPECT_EQ(mesh.mass.size(), 6);
    EXPECT_THAT(mapping, testing::ElementsAre(-1, -1, 0, 1, 2, 3, 4, 5));

    // Remove the last triangle
    ASSERT_EQ(mesh.triangles.size(), 1);
    mesh.triangles.erase(mesh.triangles.begin());
    ASSERT_EQ(mesh.triangles.size(), 0);

    mapping = remove_orphaned_vertices(mesh);
    EXPECT_EQ(mesh.vertices.size(), 0);
    EXPECT_EQ(mesh.fixed.size(), 0);
    EXPECT_EQ(mesh.mass.size(), 0);
    EXPECT_THAT(mapping, testing::ElementsAre(-1, -1, -1, -1, -1, -1));
}
}  // namespace common
