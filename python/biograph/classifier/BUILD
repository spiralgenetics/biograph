load("//tools:spiral_python.bzl", "spiral_py_library", "spiral_py_test")

spiral_py_library(
    name = "classifier",
    srcs = [
        "__init__.py",
	"build_classifier.py",
        "qual_classifier.py",
	"qual_classifier_PP.py",
	"gt_classifier_PP.py"
    ],
    lint_deps = [
        "//python/biograph:all_modules",
        "//python/biograph/classifier",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//python/biograph",
    ],
)

# spiral_py_test(
# 	name = "classifier_test",
# 	srcs = ["classifier_test.py"],
# 	deps = ["//python/biograph:all_modules"],
# )
