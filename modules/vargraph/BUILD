cc_library(
    name = "vargraph",
    srcs = ["vargraph.cpp"],
    hdrs = ["vargraph.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_test(
    name = "vargraph_test",
    srcs = ["vargraph_test.cpp"],
    deps = [
        ":vargraph",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_binary(
    name = "vargraph_perf",
    srcs = ["vargraph_perf.cpp"],
    deps = [
        ":vargraph",
    ],
)
