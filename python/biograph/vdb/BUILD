load("//tools:spiral_python.bzl", "spiral_py_library", "spiral_py_test")

spiral_py_library(
    name = "vdb",
    srcs = [
        "__init__.py",
        "anno_cmd.py",
        "athena.py",
        "cache.py",
        "filter.py",
        "parquet.py",
        "parquet_anno.py",
        "parquet_variants.py",
        "query_cmd.py",
        "sql_tables.py",
        "study_cmd.py",
        "vcf_cmd.py",
    ],
    lint_deps = ["//python/biograph:all_modules"],
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
        "//python/biograph:capi",
    ],
)

spiral_py_test(
    name = "parquet_test",
    srcs = ["parquet_test.py"],
    data = ["//golden/ftest/vdb"],
    deps = [
        "//python/biograph",
        "//python/biograph:all_modules"
    ],
)
