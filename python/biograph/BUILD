load("//tools:spiral_python.bzl", "spiral_par_binary", "spiral_py_library", "spiral_py_lint", "spiral_py_multiversion_module", "spiral_py_sdist", "spiral_py_test")

exports_files(["requirements.txt"])

cc_library(
    name = "wrappers",
    srcs = [
        "biograph.cpp",
        "biograph_metadata.cpp",
        "dna_sequence.cpp",
        "log.cpp",
        "readmap.cpp",
        "reference.cpp",
        "seqset.cpp",
        "version.cpp",
    ],
    hdrs = [
        "biograph.h",
        "biograph_metadata.h",
        "log.h",
        "readmap.h",
        "reference.h",
        "seqset.h",
    ],
    visibility = ["//python/biograph:__subpackages__"],
    deps = [
        "//base",
        "//external:pybind11",
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io",
        "//python:common",
        "//tools:libpython",
        "//vendor/htslib",
    ],
)

spiral_py_multiversion_module(
    name = "_capi",
    module = "_capi.so",
    visibility = [":__subpackages__"],
)

py_library(
    name = "capi",
    data = [":_capi"],
    visibility = ["//visibility:public"],
)

load("@rules_python//experimental/python:wheel.bzl", "py_package", "py_wheel")

cc_binary(
    name = "_capi.so",
    srcs = [
        "module.cpp",
    ],
    linkshared = 1,
    stamp = 1,
    tags = ["noasan"],
    visibility = ["//modules/bindings/python:__pkg__"],
    deps = [
        ":wrappers",
        "//external:pybind11",
        "//python/biograph/internal:wrappers",
        "//python/biograph/variants:wrappers",
        "//tools:libpython",
    ],
)

load(
    "@biograph_requirements//:requirements.bzl",
    _biograph_requirements = "all_requirements",
)

spiral_py_library(
    name = "biograph",
    srcs = [
        "__init__.py",
        "__main__.py",
        "utils.py",
    ],
    lint_deps = ["//python/biograph:all_modules"],
    visibility = ["//visibility:public"],
    deps = [":capi"] + _biograph_requirements,
)

spiral_py_library(
    name = "all_modules",
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
        "//python/biograph/classifier",
        "//python/biograph/coverage",
        "//python/biograph/internal",
        "//python/biograph/tools",
        "//python/biograph/variants",
        "//python/biograph/vdb",
    ],
)

load("//:versions.bzl", "BIOGRAPH_VERSION")

SDK_PKG = "BioGraph-" + BIOGRAPH_VERSION + ".tgz"

spiral_py_sdist(
    name = "package",
    out = SDK_PKG,
    bins = ["//modules/biograph:bgbinary"],
    packages = [
        "biograph",
        "biograph.classifier",
        "biograph.coverage",
        "biograph.internal",
        "biograph.tools",
        "biograph.variants",
        "biograph.vdb",
    ],
    readme = "//modules/bindings/python:README.txt",
    requirements = "requirements.txt",
    visibility = ["//visibility:public"],
    deps = [":all_modules"],
)

spiral_py_test(
    name = "biograph_test",
    srcs = ["biograph_test.py"],
    data = ["//golden:e_coli_bg"],
    deps = [":all_modules"],
)

spiral_py_test(
    name = "seqset_test",
    srcs = ["seqset_test.py"],
    data = ["//golden:e_coli_bg"],
    deps = [":all_modules"],
)

spiral_py_test(
    name = "readmap_test",
    srcs = ["readmap_test.py"],
    data = ["//golden:e_coli_bg"],
    deps = [":all_modules"],
)

#spiral_py_test(
#name = "log_test",
#srcs = ["log_test.py"],
#deps = [":all_modules"],
#)

spiral_py_test(
    name = "dna_test",
    srcs = ["dna_test.py"],
    deps = [":all_modules"],
)

spiral_py_test(
    name = "doctest_test",
    srcs = ["doctest_test.py"],
    data = ["//datasets"],
    tags = ["no-github-ci"],
    deps = [":all_modules"],
)
