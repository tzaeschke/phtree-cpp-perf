load("@rules_shell//shell:sh_binary.bzl", "sh_binary")

package(default_visibility = ["//visibility:private"])
load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "tinspin",
    testonly = True,
    hdrs = [
        "kdtree.h",
        "min_max_helpers.h",
        "min_max_vector_heap.h",
        "min_max_tree_heap.h",
        "quadtree_point.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//src/util:ph-util",
        "@phtree",
    ],
    alwayslink = 1,
)
