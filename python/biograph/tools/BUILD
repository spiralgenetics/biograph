load("//tools:spiral_python.bzl", "spiral_par_binary", "spiral_py_library", "spiral_py_lint")
load(
    "@biograph_requirements//:requirements.bzl",
    _biograph_requirements = "all_requirements",
)

spiral_py_library(
    name = "tools",
    srcs = [
        "__init__.py",
        "bg_qc_stats.py",
        "bgbinary_cmds.py",
        "coverage.py",
        "discover.py",
        "export_aligned.py",
        "freq_filter.py",
        "install_tests.py",
        "log.py",
        "mendelian_annotator.py",
        "merge_vcfs.py",
        "parallel_regions.py",
        "refhash.py",
        "rvg_cleaner.py",
        "vcf_sample_paste.py",
        "vcf_to_ml_table.py",
        "vdb.py",
        "visualize_tool.py",
    ],
    lint_deps = ["//python/biograph:all_modules"],
    visibility = ["//visibility:public"],
    deps = _biograph_requirements,
)
