alias(
    name = "graph_discover",
    actual = "discover",
    visibility = ["//visibility:public"],
)

cc_library(
    name = "discover",
    srcs = [
        "branch.cpp",
        "discover.cpp",
        "graph_trim_ref.cpp",
        "make_ref.cpp",
        "push_to_pair.cpp",
        "update_rc_seqset_entries.cpp",
    ],
    hdrs = [
        "branch.h",
        "discover.h",
        "graph_trim_ref.h",
        "make_ref.h",
        "push_to_pair.h",
        "update_rc_seqset_entries.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/variants:assemble",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_set",
    ],
)

cc_test(
    name = "discover_test",
    srcs = ["discover_test.cpp"],
    deps = [
        ":discover",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/test:gtest_main",
        "//modules/variants:assemble_testutil",
    ],
)

cc_test(
    name = "branch_test",
    srcs = ["branch_test.cpp"],
    deps = [
        ":discover",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/test:gtest_main",
        "//modules/variants:assemble_testutil",
    ],
)

cc_test(
    name = "push_to_pair_test",
    srcs = ["push_to_pair_test.cpp"],
    deps = [
        ":discover",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/test:gtest_main",
        "//modules/variants:assemble_testutil",
    ],
)

cc_test(
    name = "graph_trim_ref_test",
    srcs = ["graph_trim_ref_test.cpp"],
    deps = [
        ":discover",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/test:gtest_main",
        "//modules/variants:assemble_testutil",
    ],
)

cc_test(
    name = "make_ref_test",
    srcs = ["make_ref_test.cpp"],
    deps = [
        ":discover",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/test:gtest_main",
        "//modules/variants:assemble_testutil",
    ],
)
