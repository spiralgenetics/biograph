load("//tools:spiral_python.bzl", "spiral_py_library")

cc_library(
    name = "wrappers",
    srcs = [
        "anchor.cpp",
        "assemble.cpp",
        "make_vars.cpp",
        "module.cpp",
        "vargraph.cpp",
    ],
    hdrs = [
        "anchor.h",
        "assemble.h",
        "module.h",
        "vargraph.h",
    ],
    visibility = ["//python/biograph:__subpackages__"],
    deps = [
        "//external:pybind11",
        "//modules/vargraph",
        "//python:common",
        "//python/biograph/variants:wrappers",
        "//tools:libpython",
    ],
)

load(
    "@tool_requirements//:requirements.bzl",
    _tool_requirement = "requirement",
)

spiral_py_library(
    name = "internal",
    srcs = [
        "__init__.py",
        "GenomeGraph.py",
        "kPCMP.py",
        "refhashes.py",
        "validate_tool.py",
        "vPCMP.py",
    ],
    lint_deps = ["//python/biograph:all_modules"],
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
        "//python/biograph:capi",
    ],
)
