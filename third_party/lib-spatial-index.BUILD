load("@rules_cc//cc:defs.bzl", "cc_library")
licenses(["notice"])

cc_library(
    name = "lib-spatial-index",
    hdrs = glob([
#        "include/**/*.cc",
        "include/**/*.h",
    ],
     allow_empty=True),
    srcs = glob([
        "src/**/*.cc",
        "src/**/*.h",
    ],
             allow_empty=True),
    includes = ["include"],
    visibility = ["//visibility:public"],
)
