package(default_visibility = ["//visibility:private"])

#cc_library(
#    name = "pcl-kd",
#    testonly = True,
#    srcs = [
#        "src/kdtree_flann.cpp",
#    ],
#    hdrs = [
#        "include/pcl/kdtree/impl/io.hpp",
#        "include/pcl/kdtree/impl/kdtree_flann.hpp",
#        "include/pcl/kdtree/io.h",
#        "include/pcl/kdtree/kdtree.h",
#        "include/pcl/kdtree/kdtree_flann.h",
#    ],
#    visibility = [
#        "//visibility:public",
#    ],
#    deps = [
#    ],
#    alwayslink = 1,
#)

#cc_library(
#    name = "pcl-oct",
#    testonly = True,
#    srcs = [
#        "src/octree_inst.cpp",
#    ],
#    hdrs = [
#        "include/pcl/octree/octree.h",
#    ],
#    visibility = [
#        "//visibility:public",
#    ],
#    deps = [
#    ],
#    alwayslink = 1,
#)

cc_library(
    name = "ph-pcl",
    testonly = True,
    hdrs = [
        "ph-kdtree.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        #        "//src/pcl/include:pcl-kd",
        #        "//src/pcl/include:pcl-oct",
        "@pcl//:kdtree",
        "@phtree",
    ],
    alwayslink = 1,
)
