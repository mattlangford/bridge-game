cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [
        "//builder",
        "//engine:simulate",
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
        "//engine:context",
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "test",
    srcs = glob(["common/test/*.cc"]),
    deps = [
        ":common",
        "//builder:builder",
        "@gtest//:gtest",
        "@gtest//:gtest_main"
    ]
)
