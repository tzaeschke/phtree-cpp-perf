package(default_visibility = ["//visibility:private"])

cc_library(
    name = "benchmark",
    testonly = True,
    hdrs = [
        "benchmark_util.h",
        "logging.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
    alwayslink = 1,
)

cc_binary(
    name = "count_mm_d_benchmark",
    testonly = True,
    srcs = [
        "count_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "erase_benchmark",
    testonly = True,
    srcs = [
        "erase_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "erase_d_benchmark",
    testonly = True,
    srcs = [
        "erase_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "extent_benchmark",
    testonly = True,
    srcs = [
        "extent_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "extent_benchmark_weird",
    testonly = True,
    srcs = [
        "extent_benchmark_weird.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "find_benchmark",
    testonly = True,
    srcs = [
        "find_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_benchmark",
    testonly = True,
    srcs = [
        "insert_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_d_benchmark",
    testonly = True,
    srcs = [
        "insert_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/mcxme:ph-mcxme",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_mm_d_benchmark",
    testonly = True,
    srcs = [
        "insert_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/bb-tree:ph-bb-tree",
        "//src/boost:ph-boost-mm",
        "//src/flann:ph-flann",
        "//src/lsi:ph-lsi-mm",
        "//src/mcxme:ph-mcxme",
        "//src/tinspin",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_box_d_benchmark",
    testonly = True,
    srcs = [
        "insert_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/mcxme:ph-mcxme",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_mm_box_d_benchmark",
    testonly = True,
    srcs = [
        "insert_mm_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/boost:ph-boost-mm",
        "//src/lsi:ph-lsi-mm",
        "//src/mcxme:ph-mcxme",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "knn_d_benchmark",
    testonly = True,
    srcs = [
        "knn_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "knn_mm_d_benchmark",
    testonly = True,
    srcs = [
        "knn_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/boost:ph-boost-mm",
        "//src/flann:ph-flann",
        "//src/tinspin",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_benchmark",
    testonly = True,
    srcs = [
        "query_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_box_d_benchmark",
    testonly = True,
    srcs = [
        "query_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/mcxme:ph-mcxme",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_d_benchmark",
    testonly = True,
    srcs = [
        "query_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_d_benchmark",
    testonly = True,
    srcs = [
        "query_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/bb-tree:ph-bb-tree",
        "//src/boost:ph-boost-mm",
        "//src/flann:ph-flann",
        "//src/lsi:ph-lsi-mm",
        "//src/mcxme:ph-mcxme",
        "//src/robin_hood",
        "//src/tinspin",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_d_web_benchmark",
    testonly = True,
    srcs = [
        "query_mm_d_web_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/bb-tree:ph-bb-tree",
        "//src/boost:ph-boost-mm",
        "//src/flann:ph-flann",
        "//src/lsi:ph-lsi-mm",
        "//src/mcxme:ph-mcxme",
        "//src/robin_hood",
        "//src/tinspin",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_d_filter_benchmark",
    testonly = True,
    srcs = [
        "query_mm_d_filter_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_box_d_benchmark",
    testonly = True,
    srcs = [
        "query_mm_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/boost:ph-boost-mm",
        "//src/lsi:ph-lsi-mm",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_box_d_web_benchmark",
    testonly = True,
    srcs = [
        "query_mm_box_d_web_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/boost:ph-boost-mm",
        "//src/lsi:ph-lsi-mm",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_d_benchmark",
    testonly = True,
    srcs = [
        "update_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/mcxme:ph-mcxme",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_mm_d_benchmark",
    testonly = True,
    srcs = [
        "update_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "//src/boost:ph-boost-mm",
        "//src/flann:ph-flann",
        "//src/lsi:ph-lsi-mm",
        "//src/tinspin",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_mm_box_d_benchmark",
    testonly = True,
    srcs = [
        "update_mm_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_box_d_benchmark",
    testonly = True,
    srcs = [
        "update_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "hd_insert_d_benchmark",
    testonly = True,
    srcs = [
        "hd_insert_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "hd_erase_d_benchmark",
    testonly = True,
    srcs = [
        "hd_erase_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "hd_query_d_benchmark",
    testonly = True,
    srcs = [
        "hd_query_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)

cc_binary(
    name = "hd_knn_d_benchmark",
    testonly = True,
    srcs = [
        "hd_knn_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        ":benchmark",
        "@gbenchmark//:benchmark",
        "@phtree",
        "@spdlog",
    ],
)
