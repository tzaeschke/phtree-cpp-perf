package(default_visibility = ["//visibility:private"])

cc_library(
    name = "bb-tree",
    testonly = True,
    srcs = [
        "BBTree.cpp",
        "BBTreeBucket.cpp",
    ],
    hdrs = [
        "BBTree.h",
        "BBTreeBucket.h",
        "ctpl_stl.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
    ],
    alwayslink = 1,
)

cc_library(
    name = "ph-bb-tree",
    testonly = True,
    hdrs = [
        "bb_multimap.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":bb-tree",
        "@phtree",
    ],
    alwayslink = 1,
)

cc_binary(
    name = "bb-tree-main",
    testonly = True,
    srcs = [
        "main.cpp",
    ],
    linkstatic = True,
    deps = [
        ":bb-tree",
    ],
)
