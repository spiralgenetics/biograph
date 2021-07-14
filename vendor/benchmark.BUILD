cc_library(
    name = "benchmark",
    srcs = glob(
        ["src/*.cc"],
        exclude = [
            "src/re_posix.cc",
            "src/gnuregex.cc",
        ],
    ),
    hdrs = glob(
        [
            "src/*.h",
            "include/benchmark/*.h",
        ],
        exclude = [
            "src/re_posix.h",
            "src/gnuregex.h",
        ],
    ),
    copts = [
        "-DHAVE_POSIX_REGEX",
        "-DNDEBUG",
    ],
    includes = [
        "include",
    ],
    visibility = ["//visibility:public"],
)
