package(default_visibility = ["//visibility:private"])

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
