cc_library(
    name = "pipeline",
    srcs = [
        "build_reference_task.cpp",
        "client_handler.cpp",
        "dataset_meta.cpp",
        "dataset_path.cpp",
        "load_reference_dataset_meta.cpp",
        "multiple_output_step.cpp",
        "ottoman.cpp",
        "primitives.cpp",
        "steptype.cpp",
        "version.cpp",
    ],
    hdrs = [
        "all_steps.h",
        "build_reference_task.h",
        "client_handler.h",
        "dataset_meta.h",
        "dataset_path.h",
        "datatype.h",
        "direntry.h",
        "find_kmer_cutoff.h",
        "load_reference_dataset_meta.h",
        "multiple_output_step.h",
        "ottoman.h",
        "paired_merger.h",
        "primitives.h",
        "restful_registry.h",
        "steptype.h",
        "version.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//base",
        "//modules/bio_format:fasta_ref_importer",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/web:web_ssl",
        "@boost//:regex",
    ],
    alwayslink = 1,
)

cc_test(
    name = "pipeline_test",
    srcs = ["pipeline_test.cpp"],
    data = [
        "//golden:uncategorized_test_data",
    ],
    deps = [
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/pipeline",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "build_ref_test",
    srcs = ["build_ref_test.cpp"],
    data = ["//datasets"],
    deps = [
        ":unittest_primitives",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "dataset_test",
    srcs = ["dataset_test.cpp"],
    deps = [
        ":unittest_primitives",
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/pipeline",
        "//modules/test:gtest_main",
        "//modules/web",
    ],
)

cc_library(
    name = "unittest_primitives",
    testonly = 1,
    srcs = ["unittest_primitives.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":pipeline",
        "@gtest//:no_main",
    ],
    alwayslink = 1,
)
