# Depend on this to add //python to the import path.
py_library(
    name = "python_import_base",
    imports = ["."],
    visibility = ["//visibility:public"],
)

# Common C utilities for spiral python modules
cc_library(
    name = "common",
    srcs = ["common.cpp"],
    hdrs = ["common.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//external:pybind11",
        "//modules/io:base",
        "@boost//:optional",
    ],
)
