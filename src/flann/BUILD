package(default_visibility = ["//visibility:private"])

cc_library(
    name = "flann-kd",
    testonly = True,
    hdrs = [
        "flann/algorithms/dist.h",
        "flann/algorithms/kdtree_index.h",
        "flann/algorithms/nn_index.h",
        "flann/config.h",
        "flann/defines.h",
        "flann/general.h",
        "flann/util/allocator.h",
        "flann/util/any.h",
        "flann/util/dynamic_bitset.h",
        "flann/util/heap.h",
        "flann/util/matrix.h",
        "flann/util/params.h",
        "flann/util/random.h",
        "flann/util/result_set.h",
        "flann/util/saving.h",
        "flann/util/serialization.h",
    ],
    #includes = ["/src/flann"],
    strip_include_prefix = "/src/flann",
    visibility = [
        "//visibility:public",
    ],
    deps = [
    ],
    alwayslink = 1,
)

cc_library(
    name = "ph-flann",
    testonly = True,
    hdrs = [
        "ph-kdtree.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "flann-kd",
        #"@system_libs//:flann",
        "@phtree",
    ],
    alwayslink = 1,
)
