# We want to test with the following distributions:
#
# Centos 7.1, since this is what CIDRA uses
#
# Others we think are popular:
#
# Centos 8
# ubuntu >= 16.04

load("//:versions.bzl", "BIOGRAPH_VERSION")

# Change e.g. "1.2.3"-RC0-dev to just "1.2.3"
BIOGRAPH_VERSION_BASE = BIOGRAPH_VERSION.split("-")[0]

LAMBDA_TOY_DATA = "lambdaToyData-" + BIOGRAPH_VERSION_BASE + ".tgz"

MODEL = "biograph_model-" + BIOGRAPH_VERSION_BASE + ".ml"

# Repos available to build biograph
BUILD_IMAGES = [
    "manylinux2014-v2",
    # "manylinux2010",
]

# Images to run install tests on.
TEST_IMAGES = [
    # centos7.  All use python 3.6.8.
    "centos70",
    "centos71",
    "centos72",
    "centos73",
    "centos74",
    "centos75",
    "centos76",
    "centos77",
    "centos78",
    # centos7.  Python 3.6.8.
    "centos8",
    # ubuntu 14.04 (trusty), from deadsnakes.
    "ubuntu1404-py36",
    # ubuntu 16.04 (xenial), multiple python versions from deadsnakes
    "ubuntu1604-py36",
    "ubuntu1604-py37",
    "ubuntu1604-py38",
    # ubuntu 18.04 (bionic), system python 3.6.9
    "ubuntu1804",
    # ubuntu 20.04, system python 3.8.2
    "ubuntu2004",
]

load("docker.bzl", "centos7_tmpl", "docker_image_builder", "docker_pip_install_tests", "ubuntu_deadsnake_tmpl", "ubuntu_tmpl")

centos7_tmpl("centos70", "7.0.1406")

centos7_tmpl("centos71", "7.1.1503")

centos7_tmpl("centos72", "7.2.1511")

centos7_tmpl("centos73", "7.3.1611")

centos7_tmpl("centos74", "7.4.1708")

centos7_tmpl("centos75", "7.5.1804")

centos7_tmpl("centos76", "7.6.1810")

centos7_tmpl("centos77", "7.7.1908")

centos7_tmpl("centos78", "7.8.2003")

ubuntu_deadsnake_tmpl("ubuntu1404-py36", "3.6", "14.04")

# Deadsnakes ubuntu 14.04 (trusty) support:
#   python3.7 has no SSL module, so we can't download dependencies
#      (see https://github.com/deadsnakes/issues/issues/63)
#   python3.8 is not present at all.

ubuntu_deadsnake_tmpl("ubuntu1604-py36", "3.6", "16.04")

ubuntu_deadsnake_tmpl("ubuntu1604-py37", "3.7", "16.04")

ubuntu_deadsnake_tmpl("ubuntu1604-py38", "3.8", "16.04")

ubuntu_tmpl(
    "ubuntu1804",
    tag = "18.04",
)

ubuntu_tmpl(
    "ubuntu2004",
    tag = "20.04",
)

filegroup(
    name = "dockerfiles",
    srcs = [
        ":build-dockerfiles",
        ":test-dockerfiles",
    ],
)

filegroup(
    name = "build-dockerfiles",
    srcs = ["Dockerfile-" + image + "-build" for image in BUILD_IMAGES],
)

filegroup(
    name = "test-dockerfiles",
    srcs = ["Dockerfile-" + image for image in TEST_IMAGES],
)

docker_image_builder(
    name = "build",
    srcs = ["build.sh"],
    build_images = BUILD_IMAGES,
    data = [
        ":dockerfiles",
    ],
    test_images = TEST_IMAGES,
)

docker_pip_install_tests(
    images = TEST_IMAGES,
    lambda_toy_data = LAMBDA_TOY_DATA,
    model = MODEL,
    pip_package = "//python/biograph:package",
    repo = "$(DOCKER_REPO_BASE)/spiral/testbase",
    tags = ["manual"],
)

test_suite(
    name = "install_tests",
    tags = ["manual"],
    tests = [":" + image + "_install_test" for image in TEST_IMAGES],
    visibility = ["//visibility:public"],
)
