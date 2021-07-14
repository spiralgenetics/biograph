load("//tools:spiral_python.bzl", "spiral_py_library", "spiral_py_test")

spiral_py_library(
    name = "coverage",
    srcs = [
        "__init__.py",
        "anno_obj.py",
        "coverage_annotation.py",
        "genotype_annotation.py",
        "phase_conflict_resolver.py",
        "sam_output.py",
        "vcf_entry_info.py",
    ],
    lint_deps = [
        "//python/biograph:all_modules",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
        "//python/biograph/variants",
    ],
)

spiral_py_test(
    name = "covanno_test",
    srcs = ["covanno_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "gtanno_test",
    srcs = ["gtanno_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)

spiral_py_test(
    name = "sam_output_test",
    srcs = ["sam_output_test.py"],
    data = ["//datasets"],
    deps = ["//python/biograph:all_modules"],
)
