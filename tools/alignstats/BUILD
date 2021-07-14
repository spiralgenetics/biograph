cc_binary(
    name = "alignstats",
    srcs = glob([
        "*.h",
        "*.c",
    ]),
    copts = [
        "-std=gnu99",
        "-Wno-unused-result",
        "-Wno-stringop-truncation",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//vendor/htslib",
    ],
)
