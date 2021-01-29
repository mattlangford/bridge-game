cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [
        ":engine",
    ]
)

cc_library(
    name = "engine",
    hdrs = glob(["*.hh"]),
    visibility = ["//visibility:public"],
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
