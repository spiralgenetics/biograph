licenses(["notice"])

cc_library(
    name = "glog",
    srcs = [
        "src/demangle.cc",
        "src/logging.cc",
        "src/raw_logging.cc",
        "src/signalhandler.cc",
        "src/symbolize.cc",
        "src/utilities.cc",
        "src/vlog_is_on.cc",
        # Private headers
        "src/glog-config.h",
        "src/base/commandlineflags.h",
        "src/base/googleinit.h",
        "src/base/mutex.h",
        "src/glog/log_severity.h",
        "src/stacktrace.h",
        "src/stacktrace_generic-inl.h",
        "src/stacktrace_libunwind-inl.h",
        "src/stacktrace_powerpc-inl.h",
        "src/stacktrace_x86-inl.h",
        "src/stacktrace_x86_64-inl.h",
    ],
    hdrs = [
        "src/demangle.h",
        "src/glog/logging.h",
        "src/glog/raw_logging.h",
        "src/glog/stl_logging.h",
        "src/glog/vlog_is_on.h",
        "src/mock-log.h",
        "src/symbolize.h",
        "src/utilities.h",
    ],
    copts = [
        "-DNO_FRAME_POINTER",
        "-DHAVE_CONFIG_H",
        "-Isrc",
        "-Wno-sign-compare",
    ],
    includes = [
        "src",
    ],
    linkopts = [
        "-pthread",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//external:gflags",
        "//external:libunwind",
    ],
)

cc_library(
    name = "testlib",
    testonly = 1,
    hdrs = [
        "src/config_for_unittests.h",
        "src/googletest.h",
    ],
    includes = [
        "glog_build/k8/src",
    ],
)

cc_test(
    name = "demangle_unittest",
    srcs = [
        "src/demangle_unittest.cc",
    ],
    tags = ["manual"],
    deps = [
        ":glog",
        ":testlib",
    ],
)

cc_test(
    name = "logging_unittest",
    srcs = [
        "src/logging_unittest.cc",
    ],
    tags = ["manual"],
    deps = [
        ":glog",
        ":testlib",
    ],
)

cc_test(
    name = "symbolize_unittest",
    srcs = [
        "src/symbolize_unittest.cc",
    ],
    tags = ["manual"],
    deps = [
        ":glog",
        ":testlib",
    ],
)
