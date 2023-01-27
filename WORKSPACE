# Bazel bootstrapping

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive", "http_file")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

# NOTE: We make third_party/ its own bazel workspace because it allows to run `bazel build ...` without
# having all targets defined in third-party BUILD files in that directory buildable.
local_repository(
    name = "third_party",
    path = "third_party",
)

# External PH-Tree dependencies

http_archive(
    name = "spdlog",
    build_file = "@third_party//spdlog:BUILD",
    sha256 = "b38e0bbef7faac2b82fed550a0c19b0d4e7f6737d5321d4fd8f216b80f8aee8a",
    strip_prefix = "spdlog-1.5.0",
    url = "https://github.com/gabime/spdlog/archive/v1.5.0.tar.gz",
)

#https://github.com/libspatialindex/libspatialindex/archive/refs/tags/1.9.3.tar.gz
#https://github.com/libspatialindex/libspatialindex/releases/download/1.9.3/spatialindex-src-1.9.3.tar.gz
http_archive(
    name = "lib-spatial-index",
    build_file = "@third_party//lib-spatial-index:BUILD",
    sha256 = "7b44340a3edc55c11abfc453bb60f148b29f569cef9e1148583e76132e9c7379",
    strip_prefix = "libspatialindex-1.9.3",
    url = "https://github.com/libspatialindex/libspatialindex/archive/refs/tags/1.9.3.tar.gz",
)

http_archive(
    name = "gbenchmark",
    sha256 = "6132883bc8c9b0df5375b16ab520fac1a85dc9e4cf5be59480448ece74b278d4",
    strip_prefix = "benchmark-1.6.1",
    url = "https://github.com/google/benchmark/archive/v1.6.1.tar.gz",
)

http_archive(
    name = "gtest",
    sha256 = "b4870bf121ff7795ba20d20bcdd8627b8e088f2d1dab299a031c1034eddc93d5",
    strip_prefix = "googletest-release-1.11.0",
    url = "https://github.com/google/googletest/archive/release-1.11.0.tar.gz",
)

http_archive(
    name = "absl",
    sha256 = "dcf71b9cba8dc0ca9940c4b316a0c796be8fab42b070bb6b7cab62b48f0e66c4",
    strip_prefix = "abseil-cpp-20211102.0",
    url = "https://github.com/abseil/abseil-cpp/archive/20211102.0.tar.gz",
)

#http_archive(
#    name = "phtree",
#    strip_prefix = "phtree-cpp-1.3.0",
#    url = "https://github.com/tzaeschke/phtree-cpp/tree/fix/75-enable-cmake-import-of-phtree",
#)

#git_repository(
#    name = "phtree",
#    branch = "main",
##    commit = "8822dbd367eee7e3904d824b780b99009a4a9915",
#    remote = "https://github.com/tzaeschke/phtree-cpp",
#)

local_repository(
    name = "phtree",
    path = "/home/franky/work/phtree-cpp-2",
)

# Development environment tooling

BUILDIFIER_VERSION = "0.29.0"

http_file(
    name = "buildifier_linux",
    executable = True,
    sha256 = "4c985c883eafdde9c0e8cf3c8595b8bfdf32e77571c369bf8ddae83b042028d6",
    urls = ["https://github.com/bazelbuild/buildtools/releases/download/{version}/buildifier".format(version = BUILDIFIER_VERSION)],
)

http_file(
    name = "buildifier_macos",
    executable = True,
    sha256 = "9b108decaa9a624fbac65285e529994088c5d15fecc1a30866afc03a48619245",
    urls = ["https://github.com/bazelbuild/buildtools/releases/download/{version}/buildifier.mac".format(version = BUILDIFIER_VERSION)],
)

http_file(
    name = "buildifier_windows",
    executable = True,
    sha256 = "dc5d6ed5e3e0dbe9955f7606939c627af5a2be7f9bdd8814e77a22109164394f",
    urls = ["https://github.com/bazelbuild/buildtools/releases/download/{version}/buildifier.exe".format(version = BUILDIFIER_VERSION)],
)

http_archive(
    name = "bazel_compilation_database",
    sha256 = "bb1b812396e2ee36a50a13b03ae6833173ce643e8a4bd50731067d0b4e5c6e86",
    strip_prefix = "bazel-compilation-database-0.3.5",
    url = "https://github.com/grailbio/bazel-compilation-database/archive/0.3.5.tar.gz",
)

#################################################3

#load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
#all_content = """filegroup(name = "all", srcs = glob(["**"]), visibility = ["//visibility:public"])"""
#
#http_archive(
#        name = "rules_foreign_cc",
#        strip_prefix = "rules_foreign_cc-main",
#        url = "https://github.com/bazelbuild/rules_foreign_cc/archive/master.zip",
#)
#
#load("@rules_foreign_cc//:workspace_definitions.bzl", "rules_foreign_cc_dependencies")
#
#rules_foreign_cc_dependencies()
#
#http_archive(
#        name = "boost",
#        build_file_content = all_content,
#        strip_prefix = "boost_1_68_0",
#        sha256 = "da3411ea45622579d419bfda66f45cd0f8c32a181d84adfa936f5688388995cf",
#        urls = ["https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.gz"],
#)

#load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

_RULES_BOOST_COMMIT = "96e9b631f104b43a53c21c87b01ac538ad6f3b48"

http_archive(
    name = "com_github_nelhage_rules_boost",
    sha256 = "5ea00abc70cdf396a23fb53201db19ebce2837d28887a08544429d27783309ed",
    strip_prefix = "rules_boost-%s" % _RULES_BOOST_COMMIT,
    urls = [
        "https://github.com/nelhage/rules_boost/archive/%s.tar.gz" % _RULES_BOOST_COMMIT,
    ],
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()

### libSpatialIndex

#new_local_repository(
#    name = "system_libs",
#    # pkg-config --variable=libdir x11
#    path = "/usr/lib/x86_64-linux-gnu",
#    build_file_content = """
#cc_library(
#    name = "lib-spatial-index",
#    srcs = ["libspatialindex.so.6"],
#    visibility = ["//visibility:public"],
#)
#""",
#)

### libPCL
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "rules_pcl",
    sha256 = "987dd22ac4637093414ff2b9291d0e9e29f3ef156a12d5471eedaa2d5beb0a93",
    strip_prefix = "rules_pcl-1.1.0",
    url = "https://github.com/kgreenek/rules_pcl/archive/v1.1.0.tar.gz",
)

load("@rules_pcl//bzl:repositories.bzl", "pcl_repositories")

pcl_repositories()

# NOTE: This must be loaded after the call to pcl_repositories().
load("@rules_pcl//bzl:init_deps.bzl", "pcl_init_deps")

pcl_init_deps()

#new_local_repository(
#    name = "system_libs",
#    build_file_content = """
#cc_library(
#    name = "flann",
#    srcs = ["libflann.so"],
#    visibility = ["//visibility:public"],
#)
#""",
#    # pkg-config --variable=libdir x11
#    path = "/usr/lib/x86_64-linux-gnu",
#)
