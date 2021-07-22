workspace(name = "spiral")

# If the following fails, you may need to upgrade your version of bazel:
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "bazel_skylib",
    urls = [
        "https://github.com/bazelbuild/bazel-skylib/releases/download/1.0.3/bazel-skylib-1.0.3.tar.gz",
        "https://mirror.bazel.build/github.com/bazelbuild/bazel-skylib/releases/download/1.0.3/bazel-skylib-1.0.3.tar.gz",
    ],
    sha256 = "1c531376ac7e5a180e0237938a2536de0c54d93f5c278634818e0efc952dd56c",
)
load("@bazel_skylib//:workspace.bzl", "bazel_skylib_workspace")
bazel_skylib_workspace()

# libpoco 1.9.0 (modified to build with bazel)
git_repository(
    name = "poco",
    commit = "0fd7e028191280b80e65db86e8933445c4d6ad41",
    remote = "https://github.com/spiralgenetics/poco.git",
    shallow_since = "1596729823 +0000",
)

git_repository(
    name = "ncbi_cxx",
    commit = "f29ffb89db5865db25f809f33ca061855a7e2c12",
    remote = "https://github.com/spiralgenetics/ncbi_cxx.git",
)

git_repository(
    name = "boost",
    commit = "822295e3a30a218b0de80a13039462853d7620de",
    remote = "https://github.com/spiralgenetics/boost.git",
    shallow_since = "1607194149 +0000",
)

load("@boost//:submodules.bzl", "add_boost_submodules")

add_boost_submodules("@boost")

http_archive(
    name = "bzip2",
    build_file = "//vendor:bzip2.BUILD",
    # bzip.org is unavailable  as of Aug 2018:
    #    sha256 = "a2848f34fcd5d6cf47def00461fcb528a0484d8edef8208d6d2e2909dc61d9cd",
    #    url = "http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz",
    #
    # TODO(nils): Can we find a more official source for bzip2 now
    # that bzip.org disappeared?  (The checksum of the uncompressed
    # .tar file is the same, but this one is compresesd with bzip2)
    sha256 = "d70a9ccd8bdf47e302d96c69fecd54925f45d9c7b966bb4ef5f56b770960afa7",
    url = "http://http.debian.net/debian/pool/main/b/bzip2/bzip2_1.0.6.orig.tar.bz2",
)

http_archive(
    name = "libdaemon",
    build_file = "//vendor:libdaemon.BUILD",
    sha256 = "fd23eb5f6f986dcc7e708307355ba3289abe03cc381fc47a80bca4a50aa6b834",
    url = "http://0pointer.de/lennart/projects/libdaemon/libdaemon-0.14.tar.gz",
)

http_archive(
    name = "gtest",
    build_file = "//vendor:gtest.BUILD",
    sha256 = "f3ed3b58511efd272eb074a3a6d6fb79d7c2e6a0e374323d1e6bcbcc1ef141bf",
    url = "https://github.com/google/googletest/archive/release-1.8.0.zip",
)

http_archive(
    name = "openssl_1_1_1f",
    build_file = "//vendor:openssl.BUILD",
    sha256 = "186c6bfe6ecfba7a5b48c47f8a1673d0f3b0e5ba2e25602dd23b629975da3f35",
    url = "https://www.openssl.org/source/openssl-1.1.1f.tar.gz",
)

bind(
    name = "openssl",
    actual = "@openssl_1_1_1f//:openssl",
)

http_archive(
    name = "yaml",
    build_file = "//vendor:yaml-cpp.BUILD",
    sha256 = "e4d8560e163c3d875fd5d9e5542b5fd5bec810febdcba61481fe5fc4e6b1fd05",
    url = "https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.2.tar.gz",
)

git_repository(
    name = "subpar",
    commit = "74529f1df2178f07d34c72b3d270355dab2a10fc",
    remote = "https://github.com/google/subpar",
    shallow_since = "1478030877 -0700",
)

bind(
    name = "cc_toolchain",
    actual = "//tools/cpp:default-toolchain",
)

git_repository(
    name = "gflags_git",
    commit = "28f50e0fed19872e0fd50dd23ce2ee8cd759338e",
    remote = "https://github.com/gflags/gflags.git",
    shallow_since = "1548439139 +0000",
)

bind(
    name = "gflags",
    actual = "@gflags_git//:gflags",
)

bind(
    name = "libunwind",
    actual = "//vendor/libunwind:libunwind-k8",
)

new_git_repository(
    name = "gperftools",
    build_file = "//vendor:gperftools.BUILD",
    commit = "7bcb416c546ddc4354aa532c78e4ea8c8d00d75f",
    remote = "https://github.com/spiralgenetics/gperftools.git",
    shallow_since = "1580846605 +0000",
)

new_git_repository(
    name = "benchmark",
    build_file = "//vendor:benchmark.BUILD",
    commit = "f5ff6d0e0d3d00cf07bb8306548b637e98b13720",
    remote = "https://github.com/google/benchmark.git",
    shallow_since = "1489458619 -0600",
)

git_repository(
    name = "rules_python",
    commit = "6135186f93d46ab8551d9fe52bac97bf0c2de1ab",
    remote = "https://github.com/bazelbuild/rules_python.git",
    shallow_since = "1613499313 +0100",
)

load("@rules_python//python:pip.bzl", "pip_install", "pip_repositories")

# Make pylint and other tools' requirements available to depend on.
#
# TODO(nils): Is there some way we can get these pip_installs to run inside e.g.
# the manylinux2014 docker image?  Otherwise we have to explicitly specify
# the python version to use on the host, and hope that the installed version
# matches the one in the docker image.
pip_install(
    name = "tool_requirements",
    requirements = "//tools:tool_requirements.txt",
    python_interpreter = "python3.6",
)

pip_install(
    name = "biograph_requirements",
    requirements = "//python/biograph:requirements.txt",
    python_interpreter = "python3.6",
)

load("//tools:python_repository.bzl", "spiral_python_repositories")

spiral_python_repositories()

new_git_repository(
    name = "pybind11_git",
    build_file = "//vendor:pybind11.BUILD",
    commit = "80d452484c5409444b0ec19383faa84bb7a4d351",
    remote = "https://github.com/pybind/pybind11.git",
    shallow_since = "1571097444 +0200",
)

bind(
    name = "pybind11",
    actual = "@pybind11_git//:pybind11",
)

# abseil-cpp
git_repository(
    name = "com_google_absl",
    commit = "c45d1c09d517e145d722e00deea9be6c8be8dd57",
    remote = "https://github.com/abseil/abseil-cpp.git",
    shallow_since = "1588965758 -0400",
)

http_archive(
    name = "rules_cc",
    sha256 = "954b7a3efc8752da957ae193a13b9133da227bdacf5ceb111f2e11264f7e8c95",
    strip_prefix = "rules_cc-9e10b8a6db775b1ecd358d8ddd3dab379a2c29a5",
    urls = ["https://github.com/bazelbuild/rules_cc/archive/9e10b8a6db775b1ecd358d8ddd3dab379a2c29a5.zip"],
)
