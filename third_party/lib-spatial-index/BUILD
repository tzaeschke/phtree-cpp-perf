load("@rules_shell//shell:sh_binary.bzl", "sh_binary")

licenses(["notice"])

cc_library(
    name = "lib-spatial-index",
    hdrs = glob([
        "include/**/*.cc",
        "include/**/*.h",
    ]),
    srcs = glob([
        "src/**/*.cc",
        "src/**/*.h",
    ]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)
