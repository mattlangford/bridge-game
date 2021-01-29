cc_binary(
    name = "main",
    srcs = ["main.cc"],
    deps = [":engine"]
)

cc_library(
    name = "engine",
    hdrs = glob(["*.hh"]),
    visibility = ["//visibility:public"],
    deps = [],
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
