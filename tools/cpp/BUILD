package(default_visibility = ["//visibility:public"])

filegroup(
    name = "empty",
    srcs = [],
)

filegroup(
    name = "compiler_files",
    srcs = ["gcc_wrapper"],
)

load("//tools/cpp:crosstool.bzl", "cc_toolchain_config")


cc_toolchain_suite(
    name = "default-toolchain",
    toolchains = {
        "local|compiler": ":cc-compiler-local",
        "k8|clang": ":cc-clang-local",
        "local|clang": ":cc-clang-local",
        "k8": ":cc-compiler-local",
    },
)

cc_toolchain_config(name="toolchain-config",
                    cpu="k8",
                    compiler="compiler")
cc_toolchain(
    toolchain_config = ":toolchain-config",
    name = "cc-compiler-local",
    all_files = ":compiler_files",
    compiler_files = ":compiler_files",
    dwp_files = ":empty",
    linker_files = ":compiler_files",
    objcopy_files = ":empty",
    strip_files = ":empty",
    supports_param_files = 1,
    toolchain_identifier = "local_gcc",
)

cc_toolchain_config(name="clang-toolchain-config",
                    cpu="k8",
                    compiler="clang")
cc_toolchain(
    toolchain_config = ":clang-toolchain-config",
    name = "cc-clang-local",
    all_files = ":empty",
    compiler_files = ":empty",
    dwp_files = ":empty",
    linker_files = ":empty",
    objcopy_files = ":empty",
    strip_files = ":empty",
    supports_param_files = 1,
    toolchain_identifier = "local_clang",
)

# TODO(nils): bazel doesn't seem to follow this link, but instead
# follows the one in @bazel_tools//tools/cpp:malloc.  Figure out how
# to set the default malloc this way.
alias(
    name = "malloc",
    actual = select({
        "//conditions:default": ":system_malloc",
        "//tools:tcmalloc_setting": "@gperftools//:tcmalloc",
    }),
)

cc_library(
    name = "system_malloc",
)
