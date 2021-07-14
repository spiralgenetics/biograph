cc_library(
    name = "discovery",
    srcs = [
        "branch.cpp",
        "path.cpp",
        "pop_search.cpp",
        "push_search.cpp",
        "rejoin.cpp",
        "seqset_range_table.cpp",
        "state.cpp",
        "view.cpp",
        "walk_ref.cpp",
    ],
    hdrs = [
        "branch.h",
        "path.h",
        "pop_search.h",
        "push_search.h",
        "rejoin.h",
        "seqset_range_table.h",
        "state.h",
        "view.h",
        "walk_ref.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/variants:assemble",
        "@boost//:icl",
    ],
)

cc_library(
    name = "discovery_testutil",
    testonly = 1,
    srcs = ["discovery_testutil.cpp"],
    hdrs = ["discovery_testutil.h"],
    deps = [
        ":discovery",
        "//modules/bio_base:seqset_testutil",
        "//modules/variants:assemble_testutil",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "walk_ref_test",
    srcs = ["walk_ref_test.cpp"],
    deps = [
        ":discovery",
        ":discovery_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "push_search_test",
    srcs = ["push_search_test.cpp"],
    deps = [
        ":discovery",
        ":discovery_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "pop_search_test",
    srcs = ["pop_search_test.cpp"],
    deps = [
        ":discovery",
        ":discovery_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "rejoin_test",
    srcs = ["rejoin_test.cpp"],
    deps = [
        ":discovery",
        ":discovery_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "state_test",
    srcs = ["state_test.cpp"],
    deps = [
        ":discovery",
        ":discovery_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)
