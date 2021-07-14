cc_library(
    name = "minizip",
    srcs = [
        "ioapi.c",
        "unzip.c",
        "zip.c",
    ],
    hdrs = [
        "ioapi.h",
        "ioapi_buf.h",
        "ioapi_mem.h",
        "unzip.h",
        "zip.h",
    ],
    # No need for AES support.
    copts = [
        "-DNOCRYPT",
        "-DNOUNCRYPT",
    ],
    linkopts = ["-lz"],
    visibility = ["//visibility:public"],
)
