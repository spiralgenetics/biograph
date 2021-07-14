cc_library(
    name = "bio_base",
    srcs = [
        "align_astar.cpp",
        "align_multigene.cpp",
        "aligned_read.cpp",
        "biograph.cpp",
        "biograph_dir.cpp",
        "bwt_file.cpp",
        "call_structural.cpp",
        "flat_ref.cpp",
        "karyotype_compat.cpp",
        "overrep.cpp",
        "pileup.cpp",
        "readmap.cpp",
        "reference.cpp",
        "reference_assembly.cpp",
        "seq_position.cpp",
        "seqset.cpp",
        "struct_var.cpp",
        "unaligned_read.cpp",
    ],
    hdrs = [
        "align_astar.h",
        "align_multigene.h",
        "aligned_read.h",
        "astar.h",
        "biograph.h",
        "biograph_dir.h",
        "bwt_file.h",
        "call_structural.h",
        "corrected_read.h",
        "coverage_record.h",
        "flat_ref.h",
        "karyotype_compat.h",
        "kmer_counter.h",
        "overrep.h",
        "pileup.h",
        "readmap.h",
        "reference.h",
        "reference_assembly.h",
        "seq_position.h",
        "seqset.h",
        "seqset_anchor.h",
        "seqset_assemble.h",
        "seqset_bitmap.h",
        "spec_headers.h",
        "struct_var.h",
        "sv_call.h",
        "unaligned_read.h",
        "var_info.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":dna",
        "//base",
        "//modules/bio_format:exporter",
        "//modules/bio_format:fasta",
        "//modules/bio_format:importer",
        "//modules/io",
        "//modules/io:parallel",
        "//modules/io:sparse_multi",
        "//modules/io:spiral_file",
        "//vendor/htslib",
        "//vendor/msgpack",
        "@boost//:optional",
        "@boost//:range",
    ],
)

cc_library(
    name = "dna",
    srcs = [
        "dna_base.cpp",
        "dna_base_set.cpp",
        "dna_multiseq.cpp",
        "dna_sequence.cpp",
    ],
    hdrs = [
        "dna_base.h",
        "dna_base_set.h",
        "dna_multiseq.h",
        "dna_sequence.h",
        "kmer.h",
    ],
    visibility = ["//visibility:public"],
    deps = ["//modules/io"],
)

cc_library(
    name = "dna_testutil",
    testonly = 1,
    srcs = ["dna_testutil.cpp"],
    hdrs = ["dna_testutil.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":dna",
        "@boost//:regex",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "astar_test",
    srcs = ["astar_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "dna_test",
    srcs = ["dna_test.cpp"],
    deps = [
        "//modules/bio_base:dna",
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "kmer_test",
    srcs = ["kmer_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "multigene_test",
    srcs = ["multigene_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "overrep_test",
    srcs = ["overrep_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "pileup_test",
    srcs = ["pileup_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "tiny_align_test",
    srcs = ["tiny_align_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/test:gtest_main",
        "@boost//:multiprecision",
    ],
)

cc_test(
    name = "call_structural_test",
    srcs = ["call_structural_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "unaligned_read_test",
    srcs = ["unaligned_read_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "flat_ref_test",
    srcs = ["flat_ref_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "reference_wrapper_test",
    srcs = ["reference_wrapper_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "dna_testutil_test",
    srcs = ["dna_testutil_test.cpp"],
    deps = [
        ":dna_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "seqset_testutil",
    testonly = 1,
    srcs = ["seqset_testutil.cpp"],
    hdrs = ["seqset_testutil.h"],
    data = glob(["testdata/fake_ref/*"]),
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:seqset_flat",
        "//modules/bio_mapred",
        "//modules/build_seqset:builder",
        "//modules/build_seqset:repo_seq",
        "//modules/mapred",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "seqset_testutil_test",
    srcs = ["seqset_testutil_test.cpp"],
    deps = [
        ":seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "reference_testutil",
    testonly = 1,
    srcs = ["reference_testutil.cpp"],
    hdrs = ["reference_testutil.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io:unittest_config",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "reference_testutil_test",
    srcs = ["reference_testutil_test.cpp"],
    deps = [
        ":dna_testutil",
        ":reference_testutil",
        "//modules/test:gtest_main",
        "@gtest//:no_main",
    ],
)

cc_binary(
    name = "seqset_benchmark",
    srcs = ["seqset_benchmark.cpp"],
    data = [
        "testdata/hiv_test.seqset",
    ],
    deps = [
        "//modules/bio_base",
        "@benchmark",
    ],
)

cc_binary(
    name = "dna_benchmark",
    testonly = 1,
    srcs = ["dna_benchmark.cpp"],
    deps = [
        ":dna_testutil",
        "//modules/bio_base:dna",
        "@benchmark",
    ],
)

cc_test(
    name = "readmap_test",
    srcs = ["readmap_test.cpp"],
    data = [
        "//datasets",
        "//golden:e_coli_bg",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:seqset_flat",
        "//modules/bio_base:seqset_mergemap",
        "//modules/bio_base:seqset_merger",
        "//modules/bio_base:seqset_testutil",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_library(
    name = "seqset_mergemap",
    srcs = [
        "make_mergemap.cpp",
        "seqset_mergemap.cpp",
    ],
    hdrs = [
        "make_mergemap.h",
        "seqset_mergemap.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":seqset_flat",
        "//modules/io:bitcount",
        "//modules/io:spiral_file",
    ],
)

cc_test(
    name = "seqset_mergemap_test",
    srcs = ["seqset_mergemap_test.cpp"],
    deps = [
        ":seqset_mergemap",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "seqset_flat",
    srcs = [
        "seqset_flat.cpp",
    ],
    hdrs = [
        "seqset_flat.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/io:membuf",
        "//modules/io:parallel",
        "//modules/io:spiral_file",
        "//modules/test:coverage",
    ],
)

cc_library(
    name = "seqset_export",
    srcs = ["seqset_export.cpp"],
    hdrs = ["seqset_export.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":seqset_flat",
        "//modules/bio_base",
        "//modules/io:parallel",
        "//modules/io:spiral_file",
    ],
)

cc_test(
    name = "seqset_flat_test",
    srcs = ["seqset_flat_test.cpp"],
    deps = [
        ":dna_testutil",
        ":seqset_flat",
        ":seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "make_mergemap_test",
    srcs = ["make_mergemap_test.cpp"],
    deps = [
        ":dna_testutil",
        ":seqset_flat",
        ":seqset_mergemap",
        ":seqset_testutil",
        "//modules/test:gtest_main",
        "//modules/test:test_coverage",
    ],
)

cc_library(
    name = "seqset_merger",
    srcs = [
        "seqset_merger.cpp",
    ],
    hdrs = [
        "seqset_merger.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":seqset_flat",
        ":seqset_mergemap",
        "//modules/bio_base",
        "//modules/io",
        "//modules/io:membuf",
        "//modules/io:parallel",
        "//modules/io:spiral_file",
        "//modules/test:coverage",
    ],
)

cc_test(
    name = "seqset_merger_test",
    srcs = ["seqset_merger_test.cpp"],
    deps = [
        ":dna_testutil",
        ":seqset_flat",
        ":seqset_mergemap",
        ":seqset_merger",
        ":seqset_testutil",
        "//modules/test:gtest_main",
        "//modules/test:test_coverage",
    ],
)

cc_library(
    name = "shannon_entropy",
    srcs = ["shannon_entropy.cpp"],
    hdrs = ["shannon_entropy.h"],
    deps = [":bio_base"],
)

cc_test(
    name = "shannon_entropy_test",
    srcs = ["shannon_entropy_test.cpp"],
    deps = [
        ":dna_testutil",
        ":shannon_entropy",
        "//modules/test:gtest_main",
    ],
)

cc_binary(
    name = "shannon_entropy_benchmark",
    testonly = 1,
    srcs = ["shannon_entropy_benchmark.cpp"],
    deps = [
        ":dna_testutil",
        ":shannon_entropy",
        "@benchmark",
    ],
)

cc_library(
    name = "fast_read_correct",
    srcs = ["fast_read_correct.cpp"],
    hdrs = ["fast_read_correct.h"],
    visibility = ["//visibility:public"],
    deps = ["//modules/bio_base"],
)

cc_test(
    name = "fast_read_correct_test",
    timeout = "long",
    srcs = ["fast_read_correct_test.cpp"],
    deps = [
        ":dna_testutil",
        ":fast_read_correct",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "seqset_export_test",
    srcs = ["seqset_export_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":seqset_testutil",
        "//modules/bio_base:seqset_export",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "sv_call_testutil",
    testonly = 1,
    hdrs = ["sv_call_testutil.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":bio_base",
        ":dna_testutil",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "sv_call_testutil_test",
    srcs = ["sv_call_testutil_test.cpp"],
    deps = [
        ":sv_call_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "refmatch",
    srcs = ["refmatch.cpp"],
    hdrs = ["refmatch.h"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
    ],
)

cc_test(
    name = "refmatch_test",
    srcs = ["refmatch_test.cpp"],
    deps = [
        ":refmatch",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:reference_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "seqset_test",
    srcs = ["seqset_test.cpp"],
    data = ["//golden:e_coli_bg"],
    deps = [
        ":dna_testutil",
        ":seqset_testutil",
        "//modules/test:gtest_main",
    ],
)
