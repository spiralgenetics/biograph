load("//tools:spiral_python.bzl", "spiral_py_library")

spiral_py_library(
    name = "utils",
    srcs = [
        "__init__.py",
        "core.py",
        "defaults.py",
        "fasta.py",
        "fastq.py",
        "fileops.py",
        "sam.py",
        "setup.py",
        "vcf.py",
    ],
    lint_deps = ["//python/functest/utils"],
    visibility = ["//python/functest:__subpackages__"],
    deps = ["//tools:py_version"],
)
