# TODO(nils): Everything is in //modules, but //base is included by
# everything, so it should be shorter?  We should revisit the
# directory layout and make sure we don't splinter it unless we decide
# we really want to.

cc_library(
    name = "base",
    srcs = [
        "base.cpp",
        "command_line.cpp",
    ],
    hdrs = [
        "base.h",
        "command_line.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//tools:build_stamp",
        "//tools:version",
        "//vendor/glog",
        "@boost//:any",
        "@boost//:circular_buffer",
        "@boost//:dynamic_bitset",
        "@boost//:filesystem",
        "@boost//:format",
        "@boost//:multi_array",
        "@boost//:program_options",
        "@boost//:random",
        "@boost//:signals2",
        "@boost//:system",
        "@boost//:thread",
    ],
)

cc_test(
    name = "command_line_test",
    srcs = ["command_line_test.cpp"],
    deps = [
        ":base",
        "//modules/test:gtest_main",
    ],
)
