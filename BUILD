exports_files(["versions.bzl"])

test_suite(
    name = "big_tests",
    tags = ["manual"],
    tests = [
        "//modules/variants:big_tests",
        "//python/functest:big_tests",
    ],
)

test_suite(
    name = "install_tests",
    tags = ["manual"],
    tests = [
        "//tools/gendocker:install_tests",
    ],
)
