cc_library(
    name = "context",
    hdrs = ["context.hh"],
    deps = [
        "//:common",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "simulate",
    srcs = ["simulate.cc"],
    hdrs = ["simulate.hh"],
    deps = [
        "//:common",
        ":context",
        "@iterate//:iterate",
        "@glfw//:glfw",
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)
