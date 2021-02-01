
cc_library(
    name = "context",
    hdrs = ["context.hh"],
    deps = [
        "//:common",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "builder",
    srcs = ["builder.cc"],
    hdrs = ["builder.hh", "mesh_builder.hh"],
    deps = [
        "//:common",
        "//:renderer",
        ":context",
    ],
    visibility = ["//visibility:public"],
)

