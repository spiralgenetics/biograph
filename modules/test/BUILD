cc_library(
    name = "test_utils",
    testonly = 1,
    srcs = [
        "build_ref.cpp",
        "fastq_test_utils.cpp",
        "test_utils.cpp",
    ],
    hdrs = [
        "build_ref.h",
        "fastq_test_utils.h",
        "local_context.h",
        "test_utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/pipeline",
        "@gtest//:no_main",
    ],
)

cc_library(
    name = "set_ops",
    testonly = 1,
    hdrs = ["set_ops.h"],
    visibility = ["//visibility:public"],
    deps = [
        "@boost//:optional",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "set_ops_test",
    srcs = ["set_ops_test.cpp"],
    deps = [
        ":gtest_main",
        ":set_ops",
    ],
)

# Version of gtest main that installs proper handlers to get
# tracebacks on check failures.
cc_library(
    name = "gtest_main",
    testonly = 1,
    srcs = ["gtest_main.cpp"],
    # Link statically to enforce that we only get one definition of main.
    linkstatic = 1,
    visibility = ["//visibility:public"],
    deps = [
        "//base",
        "@gtest//:no_main",
    ],
    alwayslink = 1,
)

cc_library(
    name = "coverage",
    srcs = ["coverage.cpp"],
    hdrs = ["coverage.h"],
    visibility = ["//visibility:public"],
    deps = ["//base"],
)

cc_library(
    name = "test_coverage",
    testonly = 1,
    hdrs = ["test_coverage.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":coverage",
        "@gtest//:no_main",
    ],
)
