cc_library(
    name =
        "kmer_counter",
    srcs = [
        "kmer_counter.cpp",
    ],
    hdrs = [
        "kmer_count_table.h",
        "kmer_counter.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//vendor/libdivide",
    ],
)

cc_test(
    name = "kmer_counter_test",
    timeout = "long",
    srcs = ["kmer_counter_test.cpp"],
    deps = [
        ":kmer_counter",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "kmer_count_table_test",
    srcs = ["kmer_count_table_test.cpp"],
    deps = [
        ":kmer_counter",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "read_importer",
    srcs = ["read_importer.cpp"],
    hdrs = ["read_importer.h"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io",
    ],
)

cc_test(
    name = "read_importer_test",
    srcs = ["read_importer_test.cpp"],
    deps = [
        ":read_importer",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "repo_seq",
    srcs = [
        "expand.cpp",
        "part_counts.cpp",
        "part_repo.cpp",
        "repo_seq.cpp",
    ],
    hdrs = [
        "expand.h",
        "part_counts.h",
        "part_repo.h",
        "repo_seq.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io",
        "//modules/io:autostats",
    ],
)

cc_test(
    name = "repo_seq_test",
    srcs = ["repo_seq_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "part_repo_test",
    srcs = ["part_repo_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "expand_test",
    srcs = ["expand_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "correct_reads",
    srcs = [
        "correct_reads.cpp",
    ],
    hdrs = [
        "correct_reads.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":repo_seq",
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
        "@boost//:range",
    ],
)

cc_test(
    name = "correct_reads_test",
    srcs = ["correct_reads_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":correct_reads",
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "builder",
    srcs = [
        "builder.cpp",
    ],
    hdrs = [
        "builder.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":repo_seq",
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
    ],
)

cc_test(
    name = "builder_test",
    srcs = ["builder_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":builder",
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/bio_base:seqset_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "overrep_test",
    srcs = ["overrep_test.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":correct_reads",
        ":repo_seq",
        "//modules/bio_base:dna_testutil",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_binary(
    name = "repo_seq_benchmark",
    srcs = ["repo_seq_benchmark.cpp"],
    data = [
        "//etc/products:unittest.json",
    ],
    deps = [
        ":repo_seq",
        "//modules/bio_base",
        "//tools:malloc_select",
        "@benchmark",
    ],
)
