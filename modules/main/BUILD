cc_library(
    name = "main",
    srcs = glob(["*.cpp"]),
    hdrs = glob(["*.h"]),
    linkopts = [
        "-lz",
        "-ldl",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//base",
        "//modules/bio_base",
        "//modules/bio_format",
        "//modules/io",
        "//modules/mapred",
        "//modules/pipeline",
        "//modules/web",
        "@poco//:net",
    ],
)
