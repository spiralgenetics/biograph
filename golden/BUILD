filegroup(
    name = "test_qseq",
    testonly = 1,
    srcs = [
        "test.qseq",
        "test.qseq.bz2",
        "test.qseq.gz",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "e_coli_10000snp_fq",
    testonly = 1,
    srcs = [
        "e_coli_10000snp.fq",
        "e_coli_10000snp.fq.multiblock.gz",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "e_coli_bg",
    testonly = 1,
    srcs = glob([
        "e_coli_10000snp.bg/**",
        "e_coli_merged.bg/**",
    ]),
    visibility = ["//visibility:public"],
)

filegroup(
    name = "err091571",
    testonly = 1,
    srcs = [
        "ERR091571.bam",
        "ERR091571.sam",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "hiv_covar_table",
    testonly = 1,
    srcs = [
        "hiv_covar_table0",
        "hiv_covar_table1",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "yeast_10000",
    testonly = 1,
    srcs = [
        "ftest/original_qual_10000.bam",
        "ftest/original_quality_test.fq",
        "ftest/yeast_10000.bam",
        "ftest/yeast_10000.fq",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "fasta_sequences",
    testonly = 1,
    srcs = [
        "sequences.fasta",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "quick_e_coli_fq",
    testonly = 1,
    srcs = [
        "quick_e_coli.fq",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "structural_variation_txt",
    testonly = 1,
    srcs = [
        "structural_variation.txt",
    ],
    visibility = ["//visibility:public"],
)

# Test data that's yet to be categorized.  Prefer to use a more specific set of test data.
filegroup(
    name = "uncategorized_test_data",
    srcs = glob(["**/*"]),
    visibility = ["//visibility:public"],
)
