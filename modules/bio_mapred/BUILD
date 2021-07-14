cc_library(
    name = "bio_mapred",
    srcs = [
        "align_kmer.cpp",
        "align_seqset.cpp",
        "compute_coverage.cpp",
        "correct_reads_mapper.cpp",
        "display_var.cpp",
        "dna_sorter.cpp",
        "expand_mapper.cpp",
        "filter_assemblies.cpp",
        "filter_reads.cpp",
        "flatten_reducer.cpp",
        "flatten_seqset.cpp",
        "kcount_reducer.cpp",
        "kmer_filter_mapper.cpp",
        "kmer_hist_mapper.cpp",
        "kmer_set.cpp",
        "kmer_snps.cpp",
        "kmerize_bf.cpp",
        "kmerize_reads_mapper.cpp",
        "kmerize_ref.cpp",
        "kmers_to_db.cpp",
        "make_bwt.cpp",
        "make_readmap.cpp",
        "mem_seqset.cpp",
        "merge_covar.cpp",
        "merge_flat_seqset.cpp",
        "pair_reducer.cpp",
        "phred64_to_33_mapper.cpp",
        "prefix_reducer.cpp",
        "prefix_sorter.cpp",
        "read_correction.cpp",
        "read_pileup_reducer.cpp",
        "read_qual_mapper.cpp",
        "rotate_mapper.cpp",
        "seq_position_sorter.cpp",
        "seqset_assembly_data.cpp",
        "seqset_assembly_data_factory.cpp",
        "sort_expand.cpp",
        "struct_var_sorter.cpp",
        "sv_call_reducer.cpp",
        "var_id_mapper.cpp",
    ],
    hdrs = [
        "align_kmer.h",
        "align_seqset.h",
        "compute_coverage.h",
        "correct_reads.h",
        "correct_reads_mapper.h",
        "display_var.h",
        "dna_sorter.h",
        "expand_mapper.h",
        "filter_assemblies.h",
        "filter_reads.h",
        "flatten_reducer.h",
        "flatten_seqset.h",
        "kcount_reducer.h",
        "kmer_filter_mapper.h",
        "kmer_hist_mapper.h",
        "kmer_set.h",
        "kmer_snps.h",
        "kmerize_bf.h",
        "kmerize_reads_mapper.h",
        "kmerize_ref.h",
        "kmers_to_db.h",
        "make_bwt.h",
        "make_readmap.h",
        "mem_seqset.h",
        "merge_covar.h",
        "merge_flat_seqset.h",
        "overlap_graph.h",
        "pair_reducer.h",
        "phred64_to_33_mapper.h",
        "prefix_flatten.h",
        "prefix_reducer.h",
        "read_correction.h",
        "read_pileup_reducer.h",
        "read_qual_mapper.h",
        "rotate_mapper.h",
        "seq_position_sorter.h",
        "seqset_assembly_data.h",
        "seqset_assembly_data_factory.h",
        "sort_expand.h",
        "struct_var_sorter.h",
        "sv_call_reducer.h",
        "var_id_mapper.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//base",
        "//datavis/kmer_quality_report",
        "//modules/bio_base:fast_read_correct",
        "//modules/bio_base:seqset_mergemap",
        "//modules/bio_format",
        "//modules/bio_format:dna_io",
        "//modules/build_seqset:kmer_counter",
        "//modules/mapred",
        "//vendor/msgpack",
    ],
    alwayslink = 1,
)

cc_test(
    name = "read_trimmer_test",
    srcs = ["read_trimmer_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_mapred",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "phred64_to_33_test",
    srcs = ["phred64_to_33_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "kmerize_ref_test",
    srcs = ["kmerize_ref_test.cpp"],
    deps = [
        "//modules/bio_mapred",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "merge_covar_test",
    srcs = ["merge_covar_test.cpp"],
    data = [
        "//golden:hiv_covar_table",
    ],
    deps = [
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "align_kmer_test",
    srcs = ["align_kmer_test.cpp"],
    data = ["//golden:uncategorized_test_data"],
    deps = [
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "kmer_set_test",
    srcs = ["kmer_set_test.cpp"],
    data = ["//golden:uncategorized_test_data"],
    deps = [
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_binary(
    name = "kmer_set_benchmark",
    testonly = 1,
    srcs = ["kmer_set_benchmark.cpp"],
    data = [
        "//etc/products:unittest.json",
    ],
    deps = [
        "//modules/bio_mapred",
        "//modules/io",
        "@benchmark",
    ],
)

cc_test(
    name = "kmer_snps_test",
    srcs = ["kmer_snps_test.cpp"],
    deps = [
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "overlap_graph_test",
    srcs = ["overlap_graph_test.cpp"],
    data = ["//datasets"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/pipeline",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "bwt_test",
    srcs = ["bwt_test.cpp"],
    data = [
        "//datasets",
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)
