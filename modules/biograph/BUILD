cc_binary(
    name = "bgbinary",
    srcs = [
        "biograph_create.cpp",
        "biograph_info.cpp",
        "biograph_merge.cpp",
        "biograph_query.cpp",
        "biograph_variants.cpp",
        "bwt_query.cpp",
        "dump_biograph_flat.cpp",
        "export_fastq.cpp",
        "main.cpp",
        "migrate_readmap.cpp",
        "ref_to_bwt.cpp",
        "ref_to_seqset.cpp",
        "upgrade_readmap.cpp",
    ],
    linkopts = select(
        {
            "//conditions:default": [],
            "//tools:release_setting": ["-static"],
        },
    ),
    stamp = 1,
    visibility = ["//visibility:public"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_base:seqset_export",
        "//modules/bio_base:seqset_mergemap",
        "//modules/bio_base:seqset_merger",
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/build_seqset:builder",
        "//modules/build_seqset:correct_reads",
        "//modules/build_seqset:read_importer",
        "//modules/build_seqset:repo_seq",
        "//modules/io",
        "//modules/main",
        "//modules/pipeline",
        "//modules/variants:assemble",
        "//modules/variants:pipeline",
        "//modules/variants:ref_map",
        "//modules/variants:tracer",
        "//modules/variants:vcf",
        "//tools:malloc_select",
        "//vendor/htslib",
        "//vendor/msgpack",
        "@boost//:date_time",
    ],
)

cc_test(
    name = "biograph_test",
    srcs = ["biograph_test.cpp"],
    data = ["//datasets"],
    tags = ["noasan"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/pipeline",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
        "@gtest//:no_main",
    ],
)

cc_test(
    name = "biograph_merge_test",
    srcs = ["biograph_merge_test.cpp"],
    data = [
        "//datasets",
        # Requires the biograph binary
        "//modules/biograph:bgbinary",
    ],
    tags = ["noasan"],
    deps = [
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/bio_mapred",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_binary(
    name = "merge_example",
    srcs = ["merge_example.cpp"],
    stamp = 1,
    deps = [
        "//modules/bio_base:seqset_mergemap",
        "//modules/bio_base:seqset_merger",
    ],
)
