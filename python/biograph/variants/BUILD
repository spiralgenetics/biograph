load("//tools:spiral_python.bzl", "spiral_py_library", "spiral_py_test")

cc_library(
    name = "wrappers",
    srcs = [
        "add_ref.cpp",
        "align_reads.cpp",
        "apply_edges.cpp",
        "apply_graph.cpp",
        "assembly.cpp",
        "dedup_cov_reads.cpp",
        "discover.cpp",
        "filter_dup_align.cpp",
        "graph_discover.cpp",
        "limit_alleles.cpp",
        "module.cpp",
        "pair_cov.cpp",
        "pair_edge_cov.cpp",
        "par_pipeline.cpp",
        "phase.cpp",
        "pipeline.cpp",
        "place_pair_cov.cpp",
        "read_cov.cpp",
        "trim_ref.cpp",
    ],
    hdrs = [
        "add_ref.h",
        "align_reads.h",
        "apply_edges.h",
        "apply_graph.h",
        "assembly.h",
        "dedup_cov_reads.h",
        "discover.h",
        "filter_dup_align.h",
        "graph_discover.h",
        "limit_alleles.h",
        "module.h",
        "pair_cov.h",
        "pair_edge_cov.h",
        "par_pipeline.h",
        "phase.h",
        "pipeline.h",
        "pipeline_common.h",
        "place_pair_cov.h",
        "read_cov.h",
        "trim_ref.h",
    ],
    visibility = ["//python/biograph:__subpackages__"],
    deps = [
        "//external:pybind11",
        "//modules/graph_discover",
        "//modules/variants:align_reads",
        "//modules/variants:assemble",
        "//modules/variants:filter_dup_align",
        "//modules/variants:limit_alleles",
        "//modules/variants:pair_edge_cov",
        "//modules/variants:phase",
        "//modules/variants:pipeline",
        "//modules/variants:tracer",
        "//python:common",
        "//python/biograph:wrappers",
        "//tools:libpython",
        "@boost//:bind",
        "@com_google_absl//absl/strings",
    ],
)

spiral_py_library(
    name = "variants",
    srcs = ["__init__.py"],
    lint_deps = ["//python/biograph:all_modules"],
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
        "//python/biograph:capi",
    ],
)

spiral_py_test(
    name = "pair_edge_cov_test",
    srcs = ["pair_edge_cov_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "align_reads_test",
    srcs = ["align_reads_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "read_cov_test",
    srcs = ["read_cov_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "phase_test",
    srcs = ["phase_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "limit_alleles_test",
    srcs = ["limit_alleles_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "assembly_test",
    srcs = ["assembly_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "trim_ref_test",
    srcs = ["trim_ref_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "dedup_cov_reads_test",
    srcs = ["dedup_cov_reads_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "discover_test",
    srcs = ["discover_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "apply_edges_test",
    srcs = ["apply_edges_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "apply_graph_test",
    srcs = ["apply_graph_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "graph_discover_test",
    srcs = ["graph_discover_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "add_ref_test",
    srcs = ["add_ref_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)
