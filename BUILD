cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [
        ":builder",
        ":engine",
        ":renderer",
    ]
)

cc_library(
    name = "renderer",
    srcs = glob(["renderer/*.cc"]),
    hdrs = glob(["renderer/*.hh"]),
    deps = [
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    defines = ["GL_SILENCE_DEPRECATION"],
)

cc_library(
    name = "builder",
    srcs = glob(["builder/*.cc"]),
    hdrs = glob(["builder/*.hh"]),
    deps = [
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    defines = ["GL_SILENCE_DEPRECATION"],
)

cc_library(
    name = "engine",
    srcs = glob(["engine/*.cc"]),
    hdrs = glob(["engine/*.hh"]),
    deps = [
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    defines = ["GL_SILENCE_DEPRECATION"],
)

cc_test(
    name = "test",
    srcs = ["test.cc"],
    deps = [
        ":engine",
        "@gtest//:gtest",
        "@gtest//:gtest_main"
    ]
)
