cc_library(
    name = "bio_format",
    srcs = [
        "assembly.cpp",
        "corrected_reads.cpp",
        "fastq.cpp",
        "kmer_count.cpp",
        "kmer_quality_report_exporter.cpp",
        "make_vars.cpp",
        "nop_porters.cpp",
        "qseq.cpp",
        "read_qual.cpp",
        "sam_type.cpp",
        "struct_var.cpp",
        "supporting_reads.cpp",
        "vcf.cpp",
    ],
    hdrs = [
        "assembly.h",
        "corrected_reads.h",
        "fastq.h",
        "kmer_count.h",
        "kmer_quality_report_exporter.h",
        "make_vars.h",
        "nop_porters.h",
        "nucleic_acid_code.h",
        "qseq.h",
        "read_qual.h",
        "sam_type.h",
        "struct_var.h",
        "supporting_reads.h",
        "svdb_exporter.h",
        "vcf.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":exporter",
        ":importer",
        "//modules/bio_base",
        "//modules/mapred:ex_im_porter_data",
        "//modules/web:urlencode",
    ],
    alwayslink = 1,
)

cc_test(
    name = "entropy_test",
    srcs = ["entropy_test.cpp"],
    deps = [
        "//modules/bio_format",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "qseq_test",
    srcs = ["qseq_test.cpp"],
    deps = [
        "//modules/bio_format",
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "ref_importer_test",
    srcs = ["ref_importer_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_format",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "fastq_test",
    srcs = ["fastq_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "fasta_test",
    srcs = ["fasta_test.cpp"],
    data = [
        "//golden:fasta_sequences",
    ],
    deps = [
        "//modules/bio_format",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "kmer_counter_test",
    srcs = ["kmer_counter_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "make_vars_test",
    srcs = ["make_vars_test.cpp"],
    data = ["//datasets"],
    deps = [
        "//modules/bio_format",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "edit_distance_test",
    srcs = ["edit_distance_test.cpp"],
    data = ["//datasets"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io:unittest_config",
        "//modules/pipeline:unittest_primitives",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_binary(
    name = "fastq_benchmark",
    srcs = ["fastq_benchmark.cpp"],
    deps = [
        "//modules/bio_format",
        "//modules/io",
        "@benchmark",
    ],
)

cc_library(
    name = "importer",
    srcs = ["importer.cpp"],
    hdrs = ["importer.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/io",
        "//modules/io:registry",
        "//modules/io:simple_metadata",
    ],
)

cc_library(
    name = "exporter",
    srcs = ["exporter.cpp"],
    hdrs = ["exporter.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/io",
        "//modules/io:registry",
    ],
)

cc_library(
    name = "dna_io",
    srcs = [
        "dna_io.cpp",
    ],
    hdrs = [
        "dna_io.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/mapred",
    ],
)

cc_library(
    name = "fasta",
    srcs = ["fasta.cpp"],
    hdrs = ["fasta.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base:dna",
        "//modules/bio_format:exporter",
        "//modules/bio_format:importer",
        "//modules/io",
    ],
)

cc_library(
    name = "fasta_ref_importer",
    srcs = ["fasta_ref_importer.cpp"],
    hdrs = ["fasta_ref_importer.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":fasta",
        "//modules/bio_base",
        "//modules/bio_format:exporter",
        "//modules/bio_format:importer",
        "//modules/io",
        "//modules/mapred",
    ],
)
