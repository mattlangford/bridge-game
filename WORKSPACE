load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")

git_repository(
    name = "gtest",
    remote = "https://github.com/google/googletest",
    # If I use branch = "v1.10.x" I get warnings that say I should use these two lines instead
    commit = "703bd9caab50b139428cea1aaff9974ebee5742e",
    shallow_since = "1570114335 -0400",
)

new_git_repository(
    name = "glfw",
    remote = "https://github.com/glfw/glfw.git",
    commit = "8d7e5cdb49a1a5247df612157ecffdd8e68923d2",
    build_file = "@//:glfw.BUILD",
    # If I don't include this it gives warnings
    shallow_since = "1551813720 +0100",
)
