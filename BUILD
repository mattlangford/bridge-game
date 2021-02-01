cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [
        "//builder",
        ":engine",
        ":renderer",
    ]
)

cc_library(
    name = "common",
    srcs = glob(["common/*.cc"]),
    hdrs = glob(["common/*.hh"]),
    deps = [
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "renderer",
    srcs = glob(["renderer/*.cc"]),
    hdrs = glob(["renderer/*.hh"]),
    deps = [
        ":common",
        "//builder:context",
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "engine",
    srcs = glob(["engine/*.cc"]),
    hdrs = glob(["engine/*.hh"]),
    deps = [
        ":common",
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
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
