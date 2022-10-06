package(default_visibility = ["//visibility:private"])

#cc_binary(
#    name = "lib-si",
#    #srcs = ["findXfonts.c"],
#    deps = ["@system_libs//:lib-spatial-index"],
#)
#


#cc_library(
#    name = "lib-si",
#    testonly = True,
#    hdrs = [
#         "libmorton/include/morton_BMI.h",
#         "libmorton/include/morton_common.h",
#         "libmorton/include/morton_LUT_generators.h",
#   ],
#    visibility = [
#        "//visibility:public",
#    ],
#    deps = [
#        "@boost//:thread",
#    ],
#    alwayslink = 1,
#)


#cc_library(
#    name = "ph-mcxme",
#    testonly = True,
#    hdrs = [
#        "mcxme_map.h",
#   ],
#    visibility = [
#        "//visibility:public",
#    ],
#    deps = [
#        ":mcxme",
#       "@phtree//:phtree",
#    ],
#    alwayslink = 1,
#)

cc_library(
    name = "ph-lsi-mm",
    testonly = True,
    hdrs = [
        "lsi_multimap.h",
   ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        #":lib-si",
 #       "@system_libs//:lib-spatial-index",
        "@lib-spatial-index",
  #      "@lib-spatial-index//:lib-spatial-index",
       "@phtree//:phtree",
    ],
    alwayslink = 1,
)



#cc_binary(
#    name = "mcxme-main",
#    testonly = True,
#    srcs = [
#        "main.cpp",
#    ],
#    linkstatic = True,
#    deps = [
#       "//src/mcxme:ph-mcxme",
#    ],
#)
