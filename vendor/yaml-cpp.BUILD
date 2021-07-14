cc_library(
    name = "yaml",
    srcs = glob(["yaml-cpp-yaml-cpp-0.6.2/src/*.cpp"]),
    hdrs = glob([
        "yaml-cpp-yaml-cpp-0.6.2/src/*.h",
        "yaml-cpp-yaml-cpp-0.6.2/include/yaml-cpp/*.h",
        "yaml-cpp-yaml-cpp-0.6.2/include/yaml-cpp/*/*.h",
        "yaml-cpp-yaml-cpp-0.6.2/include/yaml-cpp/*/*/*.h",
    ]),
    copts = [
        "-Iexternal/yaml/yaml-cpp-yaml-cpp-0.6.2/include/",
        "-Wno-maybe-uninitialized",
    ],
    includes = ["yaml-cpp-yaml-cpp-0.6.2/include/"],
    visibility = ["//visibility:public"],
)
