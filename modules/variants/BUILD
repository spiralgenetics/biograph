cc_library(
    name = "ref_map",
    srcs = ["ref_map.cpp"],
    hdrs = ["ref_map.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_test(
    name = "ref_map_test",
    srcs = ["ref_map_test.cpp"],
    deps = [
        ":ref_map",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "pair_stats",
    srcs = ["pair_stats.cpp"],
    hdrs = ["pair_stats.h"],
    deps = [
        ":ref_map",
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_test(
    name = "pair_stats_test",
    srcs = ["pair_stats_test.cpp"],
    deps = [
        ":pair_stats",
        ":ref_map",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "assemble",
    srcs = [
        "apply_block.cpp",
        "apply_edges.cpp",
        "assemble.cpp",
        "dist_set.cpp",
        "read_set.cpp",
        "scaffold.cpp",
    ],
    hdrs = [
        "apply_block.h",
        "apply_edges.h",
        "assemble.h",
        "dist_set.h",
        "read_set.h",
        "scaffold.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":ref_map",
        "//modules/bio_base",
        "//modules/io:autostats",
        "//modules/io:ref_count",
        "@boost//:core",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
    ],
)

cc_library(
    name = "tracer_pipeline_parts",
    srcs = [
        "pair_counter.cpp",
        "sort.cpp",
    ],
    hdrs = [
        "pair_counter.h",
        "sort.h",
    ],
    deps = [
        ":assemble",
        "//modules/bio_base",
    ],
)

cc_library(
    name = "pipeline_parts",
    srcs = [
        "add_ref.cpp",
        "anchor_split.cpp",
        "calc_coverage.cpp",
        "dedup.cpp",
        "genotype.cpp",
        "normalize.cpp",
        "ploid_limit.cpp",
        "rvg_exclude.cpp",
        "simple_genotype_filter.cpp",
        "split_variants.cpp",
        "trim_ref.cpp",
    ],
    hdrs = [
        "add_ref.h",
        "anchor_split.h",
        "calc_coverage.h",
        "dedup.h",
        "genotype.h",
        "normalize.h",
        "ploid_limit.h",
        "rvg_exclude.h",
        "simple_genotype_filter.h",
        "split_variants.h",
        "trim_ref.h",
    ],
    deps = [
        ":assemble",
        ":path_group",
        ":tracer_pipeline_parts",
        "@boost//:icl",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
    ],
)

cc_library(
    name = "align",
    srcs = [
        "align.cpp",
        "anchor_drop.cpp",
    ],
    hdrs = [
        "align.h",
        "anchor_drop.h",
    ],
    deps = [":assemble"],
)

cc_library(
    name = "tracer",
    srcs = [
        "pop_tracer.cpp",
        "reversable_tracer.cpp",
        "trace_ref.cpp",
        "tracer.cpp",
    ],
    hdrs = [
        "pop_tracer.h",
        "reversable_tracer.h",
        "trace_ref.h",
        "tracer.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":align",
        ":assemble",
        ":tracer_pipeline_parts",
        "//modules/bio_base",
        "//modules/variants/discovery",
    ],
)

cc_test(
    name = "tracer_test",
    srcs = ["tracer_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":tracer",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "split_variants_test",
    srcs = ["split_variants_test.cpp"],
    # Split variants isn't used right now. TODO(nils): re-enable test
    # if we want to use it, or delete split_variants if we don't.
    tags = ["manual"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "trace_ref_test",
    srcs = ["trace_ref_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":ref_map",
        ":tracer",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "scaffold_test",
    srcs = ["scaffold_test.cpp"],
    deps = [
        ":assemble",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "align_match_test",
    srcs = ["align_match_test.cpp"],
    deps = [
        ":align",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "align_test",
    srcs = ["align_test.cpp"],
    deps = [
        ":align",
        ":assemble_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "ploid_limit_test",
    srcs = ["ploid_limit_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "calc_coverage_test",
    srcs = ["calc_coverage_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "add_ref_test",
    srcs = ["add_ref_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "genotype_test",
    srcs = ["genotype_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "dedup_test",
    srcs = ["dedup_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "assemble_test",
    srcs = ["assemble_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "sort_test",
    srcs = ["sort_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "anchor_split_test",
    srcs = ["anchor_split_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "trim_ref_test",
    srcs = ["trim_ref_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "pipeline_test",
    srcs = ["pipeline_test.cpp"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":pipeline",
        ":tracer",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "assemble_testutil",
    testonly = 1,
    srcs = ["assemble_testutil.cpp"],
    hdrs = ["assemble_testutil.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        ":tracer",
        "//modules/bio_base:seqset_testutil",
        "@gtest//:no_main",
    ],
)

cc_library(
    name = "big_assemble_testutil",
    testonly = 1,
    srcs = ["big_assemble_testutil.cpp"],
    hdrs = ["big_assemble_testutil.h"],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":assembly_dot",
        ":pipeline",
        ":tracer",
        "//modules/bio_base",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "@gtest//:no_main",
    ],
)

cc_library(
    name = "pipeline",
    srcs = ["pipeline.cpp"],
    hdrs = ["pipeline.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":align",
        ":assemble",
        ":pipeline_parts",
        ":tracer",
    ],
)

cc_test(
    name = "tinyhuman_test",
    size = "large",
    timeout = "eternal",
    srcs = ["tinyhuman_test.cpp"],
    tags = [
        "manual",
    ],
    deps = [
        ":assemble",
        ":assemble_testutil",
        ":big_assemble_testutil",
        ":pipeline",
        ":tracer",
        "//modules/bio_base",
        "//modules/bio_base:dna_testutil",
        "//modules/io",
        "//modules/main",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "hg001_downsample_test",
    size = "large",
    timeout = "eternal",
    srcs = ["hg001_downsample_test.cpp"],
    tags = [
        "manual",
    ],
    deps = [
        ":big_assemble_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "hg002_test",
    size = "large",
    timeout = "eternal",
    srcs = ["hg002_test.cpp"],
    tags = [
        "manual",
    ],
    deps = [
        ":big_assemble_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "assembly_dot",
    srcs = ["assembly_dot.cpp"],
    hdrs = ["assembly_dot.h"],
    deps = [":assemble"],
)

cc_library(
    name = "vcf",
    srcs = ["ploidless_vcf_export.cpp"],
    hdrs = ["ploidless_vcf_export.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_library(
    name = "path_group",
    srcs = ["path_group.cpp"],
    hdrs = ["path_group.h"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_test(
    name = "path_group_test",
    srcs = ["path_group_test.cpp"],
    # TODO(nils): Refactor test to work with new path_group behavior
    tags = ["manual"],
    deps = [
        ":assemble_testutil",
        ":path_group",
        "//modules/bio_base:dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "limit_alleles",
    srcs = ["limit_alleles.cpp"],
    hdrs = ["limit_alleles.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        "@boost//:icl",
    ],
)

cc_test(
    name = "limit_alleles_test",
    srcs = ["limit_alleles_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":limit_alleles",
        ":pipeline_parts",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "phase",
    srcs = ["phase.cpp"],
    hdrs = ["phase.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        "//modules/io:ref_count",
    ],
)

cc_test(
    name = "phase_test",
    srcs = ["phase_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":phase",
        ":pipeline_parts",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "filter_dup_align",
    srcs = [
        "filter_dup_align.cpp",
        "place_pair_cov.cpp",
    ],
    hdrs = [
        "filter_dup_align.h",
        "place_pair_cov.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "align_reads",
    srcs = ["align_reads.cpp"],
    hdrs = ["align_reads.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
    ],
)

cc_test(
    name = "filter_dup_align_test",
    srcs = ["filter_dup_align_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":filter_dup_align",
        ":pair_edge_cov",
        ":pipeline_parts",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "place_pair_cov_test",
    srcs = ["place_pair_cov_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":filter_dup_align",
        ":pair_edge_cov",
        ":pipeline_parts",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "pair_edge_cov",
    srcs = [
        "dedup_cov_reads.cpp",
        "pair_cov.cpp",
        "pair_edge_cov.cpp",
        "read_cov.cpp",
    ],
    hdrs = [
        "dedup_cov_reads.h",
        "pair_cov.h",
        "pair_edge_cov.h",
        "read_cov.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":assemble",
        ":path_group",
        "//modules/bio_base",
        "//modules/io",
        "//modules/io:ref_count",
        "@boost//:container",
        "@boost//:dynamic_bitset",
        "@boost//:icl",
        "@boost//:iterator",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
    ],
)

cc_test(
    name = "pair_edge_cov_test",
    srcs = ["pair_edge_cov_test.cpp"],
    data = ["//datasets"],
    deps = [
        ":assemble_testutil",
        ":pair_edge_cov",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "read_cov_test",
    srcs = ["read_cov_test.cpp"],
    data = ["//datasets"],
    deps = [
        ":assemble_testutil",
        ":pair_edge_cov",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "pair_cov_test",
    srcs = ["pair_cov_test.cpp"],
    deps = [
        ":assemble_testutil",
        ":filter_dup_align",
        ":pair_edge_cov",
        ":pipeline_parts",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

test_suite(
    name = "big_tests",
    tags = ["manual"],
    tests = [
        ":hg001_downsample_test",
        ":hg002_test",
        ":tinyhuman_test",
    ],
)

cc_binary(
    name = "read_cov_benchmark",
    srcs = ["read_cov_benchmark.cpp"],
    deps = [
        ":pair_edge_cov",
        "@benchmark",
    ],
)

cc_test(
    name = "read_set_test",
    srcs = ["read_set_test.cpp"],
    deps = [
        ":assemble",
        "//modules/test:gtest_main",
        "//modules/test:set_ops",
    ],
)

cc_test(
    name = "dist_set_test",
    srcs = ["dist_set_test.cpp"],
    deps = [
        ":assemble",
        "//modules/test:gtest_main",
        "//modules/test:set_ops",
        "@com_google_absl//absl/container:btree",
    ],
)
