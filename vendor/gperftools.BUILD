package(default_visibility = ["//visibility:public"])

licenses(["notice"])

AM_CXXFLAGS = [
    "-Wall",
    "-Wno-error=unused-result",
    "-Wno-error=unused-function",
    "-Wno-error=unused-variable",
    "-Wwrite-strings",
    "-Wno-error",
    "-Wno-sign-compare",
    "-fno-builtin-malloc",
    "-fno-builtin-free",
    "-fno-builtin-realloc",
    "-fno-builtin-calloc",
    "-fno-builtin-cfree",
    "-fno-builtin-memalign",
    "-fno-builtin-posix_memalign",
    "-fno-builtin-valloc",
    "-fno-builtin-pvalloc",
    "-fno-builtin",
    "-Iexternal/gperftools/src",
    "-g",
]

genrule(
    name = "gperftools-configure",
    srcs = glob(
        ["**/*"],
        exclude = [
            "src/config.h",
            "src/gperftools/tcmalloc.h",
            "src/gperftools/tcmalloc.cc",
        ],
    ) + ["//external:libunwind"],
    outs = [
        "src/config.h",
        "src/gperftools/tcmalloc.h",
        "src/gperftools/tcmalloc.cc",
    ],
    cmd = "" +
          "gperf=\"external/gperftools/\"; " +
          "outs=($(OUTS)); " +
          "cp -r $${gperf} build/; " +
          "( cd build; ./autogen.sh; ./configure --enable-libunwind --disable-shared) > gperftools-build.log 2>&1 || (cat gperftools-build.log; exit 1); " +
          "echo 'Patching config.h to remove ALIGNED_NEW_DELETE'; " +
          "sed '/ENABLE_ALIGNED_NEW_DELETE/d' -i build/src/config.h; " +
          "cp build/src/config.h \"$${outs[0]}\"; " +
          "cp build/src/gperftools/tcmalloc.h \"$${outs[1]}\"; " +
          "echo 'Patching tcmalloc.cc to remove large alloc nag...'; " +
          "sed 's/const int64 kDefaultLargeAllocReportThreshold =.*/const int64 kDefaultLargeAllocReportThreshold = -1;/g' $${gperf}/src/tcmalloc.cc > \"$${outs[2]}\" ;" +
          "",
)

cc_library(
    name = "tcmalloc",
    srcs = glob(
        [
            "src/third_party/valgrind.h",
            "src/*.cc",
            "src/*.h",
            "src/base/*.cc",
            "src/base/*.c",
            "src/base/*.h",
        ],
        exclude = [
            "src/debugallocation.cc",
            "src/tcmalloc.cc",
        ],
    ) + [
        ":gperftools-configure",
    ],
    hdrs = [
        "src/gperftools/heap-checker.h",
        "src/gperftools/heap-profiler.h",
        "src/gperftools/malloc_extension.h",
        "src/gperftools/malloc_extension_c.h",
        "src/gperftools/malloc_hook.h",
        "src/gperftools/malloc_hook_c.h",
        "src/gperftools/nallocx.h",
        "src/gperftools/profiler.h",
        "src/gperftools/stacktrace.h",
        "src/gperftools/tcmalloc.h",
    ],
    copts = AM_CXXFLAGS,
    includes = ["src/"],
    visibility = ["//visibility:public"],
    deps = ["//external:libunwind"],
    alwayslink = 1,
)
