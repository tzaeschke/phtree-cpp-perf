package(default_visibility = ["//visibility:private"])

cc_test(
    name = "mcxme_d_test",
    timeout = "long",
    srcs = [
        "mcxme_d_test.cc",
    ],
    linkstatic = True,
    deps = [
        "//src/mcxme:ph-mcxme",
        "@phtree//:phtree",
        "@gtest//:gtest_main",
    ],
)

cc_test(
    name = "mcxme_box_d_test",
    timeout = "long",
    srcs = [
        "mcxme_box_d_test.cc",
    ],
    linkstatic = True,
    deps = [
        "//src/mcxme:ph-mcxme",
        "@phtree//:phtree",
        "@gtest//:gtest_main",
    ],
)

cc_test(
    name = "boost_rtree_multimap_box_d_test",
    timeout = "long",
    srcs = [
        "boost_rtree_multimap_box_d_test.cc",
    ],
    linkstatic = True,
    deps = [
        "//src/boost:ph-boost-mm",
        "@phtree//:phtree",
        "@gtest//:gtest_main",
    ],
)


cc_test(
    name = "lsi_rtree_multimap_box_d_test",
    timeout = "long",
    srcs = [
        "lsi_rtree_multimap_box_d_test.cc",
    ],
    linkstatic = True,
    deps = [
        "//src/lsi:ph-lsi-mm",
        "@phtree//:phtree",
        "@gtest//:gtest_main",
    ],
)

