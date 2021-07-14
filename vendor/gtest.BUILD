package(default_visibility = ["//visibility:public"])

GMOCK = "googletest-release-1.8.0/googlemock"

GTEST = "googletest-release-1.8.0/googletest"

# Tests that want a main should depend on //modules/test:gtest_main so
# they get the proper handling for check failures.

alias(
    name = "no_main",
    actual = "gmock",
)

cc_library(
    name = "gmock",
    srcs = glob(
        [
            GMOCK + "/src/*.cc",
            GMOCK + "/src/*.h",
        ],
        exclude = [
            GMOCK + "/src/gmock-all.cc",
            GMOCK + "/src/gmock_main.cc",
        ],
    ),
    hdrs = glob([
        GMOCK + "/include/**/*.h",
    ]),
    copts = ["-Iexternal/gtest/" + GMOCK],
    includes = [
        "./googletest-release-1.8.0/googlemock/include",
    ],
    deps = [":gtest"],
)

cc_library(
    name = "gmock_main",
    srcs = [GMOCK + "/src/gmock_main.cc"],
    # Link statically to enforce we only get one definition of main.
    linkstatic = 1,
    deps = [":gmock"],
    alwayslink = 1,
)

cc_library(
    name = "gtest",
    srcs = glob(
        [
            GTEST + "/src/*.cc",
            GTEST + "/src/*.h",
        ],
        exclude = [
            GTEST + "/src/gtest-all.cc",
            GTEST + "/src/gtest_main.cc",
        ],
    ),
    hdrs = glob([
        GTEST + "/include/**/*.h",
    ]),
    copts = [
        "-Iexternal/gtest/" + GTEST,
    ],
    includes = [
        GTEST + "/include",
    ],
)

cc_library(
    name = "gtest_main",
    srcs = [GTEST + "/src/gtest_main.cc"],
    # Link statically to enforce we only get one definition of main.
    linkstatic = 1,
    deps = [":gtest"],
    alwayslink = 1,
)
