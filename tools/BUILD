exports_files([
    "dpkg.bzl",
    "run_pylint.py",
    "tool_requirements.txt",
    "make_sdist.py",
    "pylintrc",
    "timedelta.py",
])

config_setting(
    name = "gperftools",
    values = {
        "define": "GPERFTOOLS=1",
    },
    visibility = ["//visibility:public"],
)

load("//tools:make_versions.bzl", "make_version_h", "make_version_py")

make_version_h(
    name = "version_h",
    out_file = "version.h",
)

cc_library(
    name = "version",
    hdrs = [":version_h"],
    visibility = ["//visibility:public"],
)

make_version_py(
    name = "version_py",
    out_file = "py_version.py",
)

py_library(
    name = "py_version",
    srcs = ["py_version.py"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "build_stamp",
    hdrs = ["build_stamp.h"],
    linkstamp = "build_stamp.cpp",
    visibility = ["//visibility:public"],
)

genrule(
    name = "malloc_select_gen",
    outs = ["malloc_select.h"],
    cmd = """
cat > $(location malloc_select.h) << EOF
#pragma once

""" + select({
        "//conditions:default": "#undef TCMALLOC",
        "//tools:tcmalloc_setting": "#define TCMALLOC 1",
    }) + """
EOF
""",
)

cc_library(
    name = "malloc_select",
    srcs = ["malloc_select.cpp"],
    hdrs = [":malloc_select_gen"],
    visibility = ["//visibility:public"],
    deps = select({
        "//conditions:default": [],
        "//tools:tcmalloc_setting": ["@gperftools//:tcmalloc"],
    }),
    alwayslink = 1,
)

# This provides a locally generated SSL key and certificate chain that
# can be used for local testing.
genrule(
    name = "ssl_test_key",
    outs = [
        "ssl_test_key.key",
        "ssl_test_key.crt",
    ],
    cmd = """
(openssl req -x509 -newkey rsa:4096 -keyout $(location ssl_test_key.key) -out $(location ssl_test_key.key) -days 10000 -nodes -subj /CN=localhost
openssl req -new -x509 -key $(location ssl_test_key.key) -out $(location ssl_test_key.crt) -days 10000 -batch -subj /CN=localhost >> x509-err) > x509-output.txt 2>&1 || (cat x509-output.txt; exit 1)
""",
    visibility = ["//visibility:public"],
)

config_setting(
    name = "tcmalloc_setting",
    define_values = {
        "TCMALLOC": "1",
    },
    visibility = ["//visibility:public"],
)

config_setting(
    name = "release_setting",
    define_values = {
        "RELEASE": "1",
    },
    visibility = ["//visibility:public"],
)

genrule(
    name = "tcmalloc_so",
    srcs = ["@gperftools//:tcmalloc"],
    outs = ["tcmalloc.so"],
    cmd = """
for f in $(locations @gperftools//:tcmalloc)
do
   case $$f in
     *.so)  cp $$f $(location tcmalloc.so) ;;
   esac
done
        """,
)

sh_binary(
    name = "heapcheck",
    srcs = ["heapcheck.sh"],
    data = [":tcmalloc_so"],
)

load("//tools:python_repository.bzl", "PYTHON_VERSIONS")
load(
    "//tools:spiral_python.bzl",
    "spiral_par_binary",
    "spiral_py_library",
    "spiral_py_lint",
    "spiral_python_version_alias",
    "spiral_python_version_spec",
)

spiral_py_library(
    name = "timedelta",
    srcs = ["timedelta.py"],
    visibility = ["//visibility:public"],
)

spiral_python_version_spec(
    name = "python_version",
    build_setting_default = "host",
    visibility = ["//visibility:public"],
)

alias(
    name = "host-libpython",
    actual = "@spiral_python3.6//:libpython",
)

py_runtime(
    name = "host-pyruntime",
    interpreter_path = "/usr/bin/python3.6",
    python_version = "PY3",
)

spiral_python_version_alias(
    name = "libpython",
    actual = "//:libpython",  # relative to @spiral_python...
    host_actual = ":host-libpython",
    visibility = ["//visibility:public"],
)

spiral_python_version_alias(
    name = "python-runtime",
    actual = "//:runtime",  # relative to @spiral_python...
    host_actual = ":host-pyruntime",
    visibility = ["//visibility:public"],
)
