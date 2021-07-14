cc_library(
    name = "msgpack",
    srcs = glob([
        "src/*.c",
        "src/*.cpp",
    ]),
    hdrs = glob([
        "src/**/*.h",
        "src/**/*.hpp",
    ]),
    copts = ["-Wno-deprecated-declarations"],
    includes = ["src"],
    visibility = ["//visibility:public"],
)
