load("//tools:spiral_python.bzl", "spiral_par_binary", "spiral_py_library", "spiral_py_test")
load(
    "@biograph_requirements//:requirements.bzl",
    _biograph_requirement = "requirement",
)
load(
    "@tool_requirements//:requirements.bzl",
    _tool_requirement = "requirement",
)
load("//:versions.bzl", "BIOGRAPH_VERSION")
load("//:versions.bzl", "BIOGRAPH_MODEL_VERSION")

# Change e.g. "1.2.3"-RC0-dev to just "1.2.3"
BIOGRAPH_VERSION_BASE = BIOGRAPH_VERSION.split("-")[0]

spiral_py_test(
    name = "basic_test",
    timeout = "short",
    srcs = ["basic_test.py"],
    data = [
        "//modules/biograph:bgbinary",
    ],
    tags = ["noasan"],
    deps = [
        "//python/functest/utils",
    ],
)

spiral_py_test(
    name = "sam_test",
    timeout = "short",
    srcs = ["sam_test.py"],
    data = ["//golden:uncategorized_test_data"],
    deps = [
        "//python/functest/utils",
        _biograph_requirement("pysam"),
    ],
)

spiral_py_test(
    name = "biograph_test",
    size = "enormous",
    srcs = [
        "biograph_test.py",
    ],
    data = [
        "//datasets",
        "//golden:uncategorized_test_data",
        "//modules/biograph:bgbinary",
    ],
    tags = [
        "no-github-ci",
        "noasan",
    ],
    deps = [
        "//python/biograph:all_modules",
        "//python/functest/utils",
        _biograph_requirement("pyvcf"),
    ],
)

spiral_py_test(
    name = "sdk_test",
    # If they aren't cached in RAM, this can take a long time to read
    # the tens of gigabytes of input files.
    size = "enormous",
    timeout = "eternal",
    srcs = [
        "sdk_test.py",
    ],
    data = [
        ":biograph_main",
        #":biograph_model",
        "//datasets",
        "//golden:uncategorized_test_data",
        "//modules/biograph:bgbinary",
    ],
    tags = [
        "no-github-ci",
        "noasan",
    ],
    deps = [
        "//python/biograph:all_modules",
        "//python/functest/utils",
    ],
)

spiral_py_test(
    name = "GenomeGraphTests",
    timeout = "eternal",
    srcs = [
        "GenomeGraphTests.py",
    ],
    tags = [
        "manual",
        "noasan",
    ],
    deps = [
        "//python/biograph",
        "//python/biograph/internal",
        "//python/functest/utils",
        _biograph_requirement("pyvcf"),
        _biograph_requirement("pytabix"),
    ],
)

spiral_py_test(
    name = "import_formats_test",
    size = "enormous",
    srcs = [
        "import_formats_test.py",
    ],
    data = [
        "//datasets",
        "//modules/biograph:bgbinary",
    ],
    deps = [
        "//python/biograph",
        "//python/biograph/internal",
        "//python/functest/utils",
    ],
)

spiral_py_test(
    name = "sv_perf_check",
    size = "enormous",
    srcs = [
        "sv_perf_check.py",
    ],
    data = [
        "sv_perf_baseline.txt",
        "//modules/biograph:bgbinary",
        "//modules/nexus:biograph_create/resources/usr/bin/bgzip",
        "//modules/nexus:biograph_create/resources/usr/bin/tabix",
        "//modules/nexus:biograph_create/resources/usr/bin/vcf-sort",
        "//python/functest/utils",
    ],
    tags = [
        "manual",
    ],
    deps = [
        "//python/functest/utils",
        _tool_requirement("truvari"),
    ],
)

# Like sv_perf_check, but for the full pipeline.
spiral_py_test(
    name = "sv_perf_check_full",
    timeout = "eternal",
    srcs = [
        "sv_perf_check_full.py",
    ],
    args = [BIOGRAPH_VERSION_BASE],
    data = [
        "sv_perf_check_full.bed",
        "sv_perf_full_baseline.txt",
        ":biograph_main",
        "//modules/biograph:bgbinary",
        "//modules/nexus:biograph_create/resources/usr/bin/bgzip",
        "//modules/nexus:biograph_create/resources/usr/bin/tabix",
        "//modules/nexus:biograph_create/resources/usr/bin/vcf-sort",
    ],
    # TODO(nils): Truvari's Levenshtein dependency doesn't seem to get installed correctly
    # when this test doesn't run locally; fix this.
    local = 1,
    tags = [
        "manual",
    ],
    deps = [
        "//python/biograph:all_modules",
        "//python/functest/utils",
        _tool_requirement("truvari"),
    ],
)

spiral_py_test(
    name = "refhash_test",
    timeout = "short",
    srcs = ["refhash_test.py"],
    data = [
        "//datasets",
        "//golden:uncategorized_test_data",
        "//golden/ftest/vdb",
    ],
    tags = ["no-github-ci"],
    deps = [
        "//python/biograph:all_modules",
        "//python/functest/utils",
    ],
)

test_suite(
    name = "big_tests",
    tags = ["manual"],
    tests = [
        ":sv_perf_check",
        # TODO(nils): Re-enable this once it successfully runs under jenkins
        # with --config=manylinux2014
        #        ":sv_perf_check_full",
    ],
)

spiral_py_test(
    name = "vdb_test",
    timeout = "eternal",
    srcs = ["vdb_test.py"],
    data = [
        ":biograph_main",
        "//golden/ftest/vdb",
    ],
    tags = [
        "manual",
        "noasan",
    ],
    deps = [
        "//python/biograph",
        "//python/biograph:all_modules",
        "//python/functest/utils",
    ],
)

spiral_py_test(
    name = "filter_test",
    srcs = [
        "filter_test.py",
    ],
    data = [
        ":biograph_main",
    ],
    tags = [
        "noasan",
    ],
    deps = [
        "//python/biograph:all_modules",
        "//python/functest/utils",
    ],
)

genrule(
    name = "biograph_model",
    outs = ["biograph_model.ml"],
    cmd = """
curl -f -o $@ http://spiral-public.s3.us-west-2.amazonaws.com/models/biograph_model-""" + BIOGRAPH_MODEL_VERSION + """.ml
""",
    tags = ["no-github-ci"],
)
