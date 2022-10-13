package(default_visibility = ["//visibility:private"])

cc_library(
    name = "mcxme",
    testonly = True,
    hdrs = [
        "Entry.h",
        "PHTree.h",
        "iterators/AHCIterator.h",
        "iterators/LHCIterator.h",
        "iterators/NodeIterator.h",
        "iterators/RangeQueryIterator.h",
        "iterators/RangeQueryStackContent.h",
        "libmorton/include/morton.h",
        "libmorton/include/morton2D.h",
        "libmorton/include/morton2D_LUTs.h",
        "libmorton/include/morton3D.h",
        "libmorton/include/morton3D_LUTs.h",
        "libmorton/include/morton_BMI.h",
        "libmorton/include/morton_LUT_generators.h",
        "libmorton/include/morton_common.h",
        "nodes/AHC.h",
        "nodes/LHC.h",
        "nodes/Node.h",
        "nodes/NodeAddressContent.h",
        "nodes/SuffixStorage.h",
        "nodes/TNode.h",
        "nodes/TSuffixStorage.h",
        "util/DeletedNodes.h",
        "util/DynamicNodeOperationsUtil.h",
        "util/EntryBuffer.h",
        "util/EntryBufferPool.h",
        "util/EntryTreeMap.h",
        "util/FileInputUtil.h",
        "util/InsertionThreadPool.h",
        "util/MultiDimBitset.h",
        "util/NodeTypeUtil.h",
        "util/PlotUtil.h",
        "util/RandUtil.h",
        "util/RangeQueryThreadPool.h",
        "util/RangeQueryUtil.h",
        "util/ResultStorage.h",
        "util/SpatialSelectionOperationsUtil.h",
        "util/TEntryBuffer.h",
        "util/compare/ParallelRangeQueryScan.h",
        "util/compare/RTreeBulkWrapper.h",
        "util/rdtsc.h",
        "visitors/AssertionVisitor.h",
        "visitors/CountNodeTypesVisitor.h",
        "visitors/PrefixSharingVisitor.h",
        "visitors/SizeVisitor.h",
        "visitors/SuffixVisitor.h",
        "visitors/Visitor.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "@boost//:thread",
    ],
    alwayslink = 1,
)

cc_library(
    name = "ph-mcxme",
    testonly = True,
    hdrs = [
        "mcxme_map.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":mcxme",
        "@phtree",
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
