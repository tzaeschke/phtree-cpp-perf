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

