cc_library(
    name = "bzlib",
    srcs = [
        "bzip2-1.0.6/blocksort.c",
        "bzip2-1.0.6/bzlib.c",
        "bzip2-1.0.6/bzlib_private.h",
        "bzip2-1.0.6/compress.c",
        "bzip2-1.0.6/crctable.c",
        "bzip2-1.0.6/decompress.c",
        "bzip2-1.0.6/huffman.c",
        "bzip2-1.0.6/randtable.c",
    ],
    hdrs = ["bzip2-1.0.6/bzlib.h"],
    copts = ["-Wno-maybe-uninitialized"],
    includes = ["bzip2-1.0.6"],
    visibility = ["//visibility:public"],
)
