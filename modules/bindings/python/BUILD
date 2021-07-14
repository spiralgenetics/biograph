load("//tools:spiral_python.bzl", "spiral_py_library", "spiral_par_binary", "spiral_py_binary")

exports_files([
    "README.txt",
])

# Dependency that adds an import path in this directory so we can
# specify relative imports in the examples.
py_library(
    name = "import_path",
    imports = ["."],
)

load("//tools:release.bzl", "spiral_release")
load("//:versions.bzl", "BIOGRAPH_VERSION")
