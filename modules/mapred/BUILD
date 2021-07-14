cc_library(
    name = "mapred",
    srcs = [
        "blob_create.cpp",
        "dual_map_task.cpp",
        "dual_mapper.cpp",
        "file_info_reader.cpp",
        "histogram_export.cpp",
        "identity_mapper.cpp",
        "identity_reducer.cpp",
        "kv_hold.cpp",
        "kv_merge.cpp",
        "kv_sort.cpp",
        "kv_summarize.cpp",
        "lexical_sorter.cpp",
        "manifest.cpp",
        "map_pipe_task.cpp",
        "map_reduce_task.cpp",
        "map_task.cpp",
        "mapper.cpp",
        "mapred_wrapper_task.cpp",
        "merge_reader.cpp",
        "metadata.cpp",
        "migrate_000_001.cpp",
        "migrate_001_002.cpp",
        "migrate_002_003.cpp",
        "migrate_003_004.cpp",
        "output_stream.cpp",
        "pair_sorter.cpp",
        "path.cpp",
        "path_file.cpp",
        "path_s3_stub.cpp",
        "pipe_mapper.cpp",
        "query.cpp",
        "reduce_task.cpp",
        "reducer.cpp",
        "resource_manager.cpp",
        "run_wrapper_task.cpp",
        "sort_task.cpp",
        "sorter.cpp",
        "splitter.cpp",
        "splitter_task.cpp",
        "sum_reducer.cpp",
        "task.cpp",
        "task_composite.cpp",
        "task_info.cpp",
        "task_mgr.cpp",
        "task_runner.cpp",
        "task_tree.cpp",
        "task_worker.cpp",
        "taskdb.cpp",
        "temp_file.cpp",
        "unix_pipeline.cpp",
        "value_count_mapper.cpp",
    ],
    hdrs = [
        "base_chunker.h",
        "blob_create.h",
        "dual_map_task.h",
        "dual_mapper.h",
        "file_info_reader.h",
        "histogram_export.h",
        "identity_mapper.h",
        "identity_reducer.h",
        "input_stream.h",
        "kv_cache.h",
        "kv_hold.h",
        "kv_merge.h",
        "kv_sort.h",
        "kv_summarize.h",
        "lexical_sorter.h",
        "manifest.h",
        "manifest_parallel.h",
        "map_pipe_task.h",
        "map_reduce_task.h",
        "map_task.h",
        "mapper.h",
        "mapred_wrapper_task.h",
        "merge_reader.h",
        "metadata.h",
        "migrate.h",
        "multi_reader.h",
        "name_generator.h",
        "output_stream.h",
        "pair_sorter.h",
        "path.h",
        "path_impl.h",
        "path_s3_stub.h",
        "pipe_mapper.h",
        "pipe_params.h",
        "query.h",
        "reduce_task.h",
        "reducer.h",
        "resource_manager.h",
        "run_wrapper_task.h",
        "sort_task.h",
        "sorter.h",
        "splitter.h",
        "splitter_task.h",
        "sum_reducer.h",
        "task.h",
        "task_attempt.h",
        "task_composite.h",
        "task_info.h",
        "task_mgr.h",
        "task_runner.h",
        "task_tree.h",
        "task_worker.h",
        "taskdb.h",
        "temp_file.h",
        "unix_pipeline.h",
        "value_count_mapper.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":ex_im_porter_data",
        "//base",
        "//modules/bio_base",
        "//modules/bio_format:exporter",
        "//modules/bio_format:importer",
        "//modules/io",
        "//modules/io:registry",
        "//modules/io:simple_metadata",
        "//modules/web",
        "@boost//:integer",
        "@boost//:regex",
        "@boost//:tokenizer",
    ],
    alwayslink = 1,
)

cc_test(
    name = "manifest_test",
    srcs = ["manifest_test.cpp"],
    data = ["//datasets"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "metadata_test",
    srcs = ["metadata_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "sum_reducer_test",
    srcs = ["sum_reducer_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "chunker_test",
    srcs = ["chunker_test.cpp"],
    deps = [
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "encoding_test",
    srcs = ["encoding_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "file_wrapper_test",
    srcs = ["file_wrapper_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "find_test",
    srcs = ["find_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "keyvalue_test",
    srcs = ["keyvalue_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "merge_sort_test",
    srcs = ["merge_sort_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "mmap_buffer_test",
    srcs = ["mmap_buffer_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "pipe_mapper_test",
    srcs = ["pipe_mapper_test.cpp"],
    copts = ["-DGTEST_HAS_TR1_TUPLE=0"],
    deps = [
        "//modules/bio_format",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "progress_test",
    srcs = ["progress_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "sort_reduce_test",
    srcs = ["sort_reduce_test.cpp"],
    deps = [
        "//modules/bio_base",
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "sort_test",
    srcs = ["sort_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "task_test",
    srcs = ["task_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "task_tree_test",
    srcs = ["task_tree_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "resource_manager_test",
    srcs = ["resource_manager_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "path_test",
    srcs = ["path_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/io:unittest_config",
        "//modules/mapred",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_library(
    name = "ex_im_porter_data",
    hdrs = ["ex_im_porter_data.h"],
    visibility = ["//visibility:public"],
    deps = ["//modules/io"],
)
