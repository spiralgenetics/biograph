package(default_visibility = ["//visibility:public"])

cc_library(
    name = "io",
    deps = [
        ":base",
        ":bitcount",
        ":json",
        ":membuf",
        ":misc",
        ":mmap_file",
        ":packed_vector",
        ":runtime_stats",
        ":spiral_file",
        ":uuid",
        "//base",
        "//external:openssl",
    ],
)

cc_library(
    name = "base",
    srcs = [
        "io.cpp",
        "log.cpp",
        "track_mem.cpp",
        "utils.cpp",
        "version.cpp",
    ],
    hdrs = [
        "io.h",
        "log.h",
        "make_unique.h",
        "progress.h",
        "stopwatch.h",
        "string_view.h",
        "track_mem.h",
        "transfer_object.h",
        "utils.h",
        "version.h",
    ],
    deps = [
        "//base",
        "//tools:build_stamp",
        "//tools:version",
        "@boost//:algorithm",
        "@boost//:container",
        "@boost//:format",
        "@boost//:lexical_cast",
        "@boost//:regex",
        "@boost//:system",
        "@boost//:uuid",
    ],
)

cc_library(
    name = "parallel",
    srcs = ["parallel.cpp"],
    hdrs = ["parallel.h"],
    deps = [":base"],
)

cc_library(
    name = "json",
    srcs = [
        "json_reader.cpp",
        "json_transfer.cpp",
        "json_value.cpp",
        "json_writer.cpp",
    ],
    hdrs = [
        "json_transfer.h",
    ],
    deps = [
        ":base",
        "//vendor/json_spirit",
        "@boost//:thread",
    ],
)

cc_library(
    name = "membuf",
    srcs = [
        "membuf.cpp",
    ],
    hdrs = ["membuf.h"],
    deps = [
        ":base",
        ":parallel",
    ],
)

cc_library(
    name = "spiral_file",
    srcs = [
        "spiral_file.cpp",
        "spiral_file_mem.cpp",
        "spiral_file_mmap.cpp",
    ],
    hdrs = [
        "spiral_file.h",
        "spiral_file_mem.h",
        "spiral_file_mmap.h",
    ],
    deps = [
        ":base",
        ":json",
        ":membuf",
        ":mmap_file",
        ":uuid",
        "//vendor/minizip",
    ],
)

cc_library(
    name = "mmap_file",
    srcs = [
        "file_wrapper.cpp",
        "mmap_buffer.cpp",
    ],
    hdrs = [
        "file_wrapper.h",
        "mmap_buffer.h",
    ],
    deps = [
        ":base",
        ":membuf",
        "@boost//:format",
    ],
)

cc_library(
    name = "bitcount",
    srcs = ["bitcount.cpp"],
    hdrs = [
        "bitcount.h",
    ],
    deps = [
        ":base",
        ":packed_vector",
        ":spiral_file",
        "//base",
        "//modules/test:coverage",
    ],
)

cc_library(
    name = "packed_vector",
    srcs = [
        "int_map_interface.cpp",
        "packed_varbit_vector.cpp",
        "packed_varint_vector.cpp",
    ],
    hdrs = [
        "int_map_interface.h",
        "packed_varbit_vector.h",
        "packed_varint_vector.h",
        "packed_vector.h",
    ],
    deps = [
        ":membuf",
        ":spiral_file",
        "//base",
    ],
)

cc_test(
    name = "packed_varbit_vector_test",
    srcs = ["packed_varbit_vector_test.cpp"],
    deps = [
        ":packed_vector",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "int_map_interface_test",
    srcs = ["int_map_interface_test.cpp"],
    deps = [
        ":packed_vector",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "uuid",
    srcs = ["uuid.cpp"],
    hdrs = ["uuid.h"],
)

cc_test(
    name = "track_mem_test",
    srcs = ["track_mem_test.cpp"],
    deps = [
        ":base",
        ":membuf",
        ":misc",
        "//modules/test:gtest_main",
        "//tools:malloc_select",
    ],
)

cc_test(
    name = "uuid_test",
    srcs = ["uuid_test.cpp"],
    deps = [
        ":uuid",
        "//modules/test:gtest_main",
        "@boost//:regex",
    ],
)

cc_library(
    name = "runtime_stats",
    srcs = ["runtime_stats.cpp"],
    hdrs = [
        "runtime_stats.h",
    ],
    deps = [
        ":base",
        ":json",
        "@boost//:filesystem",
    ],
)

# TODO(nils): Break this big library up.
cc_library(
    name = "misc",
    srcs = [
        "base64.cpp",
        "bzip.cpp",
        "color_text_buffer.cpp",
        "command.cpp",
        "config.cpp",
        "crypto.cpp",
        "dict_codec.cpp",
        "digest.cpp",
        "dynamic_codecs.cpp",
        "encoding.cpp",
        "file_io.cpp",
        "hexdump.cpp",
        "inverse_kvread.cpp",
        "keyvalue.cpp",
        "pass_thru.cpp",
        "prefix_sum.cpp",
        "progress_tracker.cpp",
        "readable_PRNG.cpp",
        "tunstall.cpp",
        "zip.cpp",
    ],
    hdrs = [
        "aggregate_map.h",
        "base64.h",
        "binary_find.h",
        "blocking_queue.h",
        "bloom_filter.h",
        "bzip.h",
        "color_text_buffer.h",
        "command.h",
        "config.h",
        "crypto.h",
        "defaults.h",
        "dict_codec.h",
        "digest.h",
        "double_pair.h",
        "dynamic_codecs.h",
        "encoding.h",
        "file_io.h",
        "hash.h",
        "hash_io.h",
        "hexdump.h",
        "int_seq.h",
        "inverse_kvread.h",
        "keyvalue.h",
        "loop_io.h",
        "mem_io.h",
        "mmap_vector.h",
        "msgpack_transfer.h",
        "pass_thru.h",
        "prefix_sum.h",
        "progress_tracker.h",
        "progress_tracker_types.h",
        "pulse.h",
        "range_coder.h",
        "readable_PRNG.h",
        "semaphore.h",
        "skip_iterator.h",
        "stats.h",
        "tab_io.h",
        "thread_barrier.h",
        "tunstall.h",
        "zip.h",
    ],
    linkopts = [
        "-lz",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":base",
        ":bitcount",
        ":json",
        ":mmap_file",
        ":packed_vector",
        "//base",
        "//external:openssl",
        "//tools:build_stamp",
        "//tools:version",
        "//vendor/json_spirit",
        "//vendor/msgpack",
        "@boost//:regex",
        "@bzip2//:bzlib",
    ],
)

cc_library(
    name = "unittest_config",
    testonly = 1,
    srcs = ["unittest_config.cpp"],
    data = ["//etc/products:unittest.json"],
    visibility = ["//visibility:public"],
    deps = [
        ":io",
        "//modules/web",
        "@gtest//:no_main",
    ],
    alwayslink = 1,
)

cc_test(
    name = "io_test",
    srcs = ["io_test.cpp"],
    data = [
        "//golden:e_coli_10000snp_fq",
    ],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "aggregate_map_test",
    srcs = ["aggregate_map_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "bitcount_test",
    srcs = ["bitcount_test.cpp"],
    deps = [
        "//modules/io:bitcount",
        "//modules/test:gtest_main",
        "//modules/test:test_coverage",
    ],
)

cc_test(
    name = "bzip_test",
    srcs = ["bzip_test.cpp"],
    data = [
        "//golden:test_qseq",
    ],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "command_test",
    srcs = ["command_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "crypto_test",
    srcs = ["crypto_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "mem_io_test",
    srcs = ["mem_io_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "prefix_sum_test",
    srcs = ["prefix_sum_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "range_coder_test",
    srcs = ["range_coder_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "readable_prng_test",
    srcs = ["readable_prng_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "serialize_test",
    srcs = ["serialize_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "tab_io_test",
    srcs = ["tab_io_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "test_config",
    srcs = [
        "test_config_test.cpp",
    ],
    data = [
        "//etc/products:unittest.json",
    ],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "tunstall_test",
    srcs = ["tunstall_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "utils_test",
    srcs = ["utils_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "version_test",
    srcs = ["version_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "zip_test",
    srcs = ["zip_test.cpp"],
    data = [
        "//golden:e_coli_10000snp_fq",
        "//golden:test_qseq",
    ],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "color_text_buffer_test",
    srcs = ["color_text_buffer_test.cpp"],
    deps = [
        "//modules/io",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "packed_vector_test",
    srcs = ["packed_vector_test.cpp"],
    deps = [
        ":packed_vector",
        "//modules/test:gtest_main",
        "//modules/test:test_utils",
    ],
)

cc_test(
    name = "bloom_test",
    srcs = ["bloom_test.cpp"],
    data = ["//golden:uncategorized_test_data"],
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

cc_test(
    name = "membuf_test",
    srcs = ["membuf_test.cpp"],
    deps = [
        "//modules/io:membuf",
        "//modules/io:mmap_file",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "spiral_file_test",
    srcs = ["spiral_file_test.cpp"],
    deps = [
        "//base",
        "//modules/io:spiral_file",
        "//modules/io:unittest_config",
        "//modules/test:gtest_main",
    ],
)

cc_test(
    name = "parallel_test",
    srcs = ["parallel_test.cpp"],
    deps = [
        ":parallel",
        "//modules/test:gtest_main",
    ],
)

cc_binary(
    name = "bitcount_benchmark",
    srcs = ["bitcount_benchmark.cpp"],
    deps = [
        "//modules/io",
        "@benchmark",
    ],
)

cc_binary(
    name = "packed_varbit_vector_benchmark",
    srcs = ["packed_varbit_vector_benchmark.cpp"],
    deps = [
        "//modules/io",
        "@benchmark",
    ],
)

cc_binary(
    name = "parallel_benchmark",
    srcs = ["parallel_benchmark.cpp"],
    deps = [
        "//modules/io",
        "@benchmark",
    ],
)

cc_library(
    name = "sparse_multi",
    srcs = ["sparse_multi.cpp"],
    hdrs = ["sparse_multi.h"],
    deps = [
        ":base",
        ":bitcount",
        ":spiral_file",
    ],
)

cc_test(
    name = "sparse_multi_test",
    srcs = ["sparse_multi_test.cpp"],
    deps = [
        ":sparse_multi",
        ":spiral_file",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "autostats",
    srcs = ["autostats.cpp"],
    hdrs = ["autostats.h"],
    visibility = ["//visibility:public"],
    deps = ["@boost//:preprocessor"],
)

cc_test(
    name = "autostats_test",
    srcs = ["autostats_test.cpp"],
    deps = [
        ":autostats",
        "//modules/test:gtest_main",
    ],
)

cc_library(
    name = "registry",
    srcs = ["registry.cpp"],
    hdrs = ["registry.h"],
    visibility = ["//visibility:public"],
    deps = [":base"],
)

cc_library(
    name = "simple_metadata",
    srcs = ["simple_metadata.cpp"],
    hdrs = ["simple_metadata.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":base",
        ":json",
    ],
)

cc_library(
    name = "ref_count",
    hdrs = ["ref_count.h"],
    visibility = ["//visibility:public"],
    deps = [":base"],
)
