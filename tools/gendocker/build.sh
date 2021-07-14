#!/bin/sh

# To push build and test images to spiral docker, run this script with:
#   bazel run //tools/gendocker:build

set -e

usage() {
    cat 1>&2 <<'USAGE'
Normal Usage:
   bazel run //tools/gendocker:build -- [<options>] [<image> ...]

<options> can be:
   --push:  Push all images to repository after building.

If 1 or more <image>s are specified, only those images are built/pushed.  If no
<image> is specified, all images are built and pushed.

When this script invoked by bazel, bazel also passes the following options:
   --docker-repo-base <repo base>: Specify repo base.
   --build-image <build image name>: Specify we have a build image with the given name.
   --test-image <test image name>: Specify that we have a test image with the given name.

USAGE
    echo "$*" 1>& 2
    exit 1
}

DO_PUSH=false
arg_with_value(){
    if [ "${argsleft}" -eq 0 ]
    then
	usage "Missing argument for ${key}"
    fi
}

while [ $# -gt 0 ]
do
    key="$1"
    shift
    value="$1"
    argsleft="$#"
    case "${key}" in
	# Supplied by user:
	--push)
	    DO_PUSH=true
	    ;;
	# Supplied by bazel:
	--docker-repo-base)
	    arg_with_value
	    shift
	    DOCKER_REPO_BASE="${value}"
	    ;;
	--test-image)
	    arg_with_value
	    shift
	    TEST_IMAGES="${TEST_IMAGES} ${value}"
	    ;;
	--build-image)
	    arg_with_value
	    shift
	    BUILD_IMAGES="${BUILD_IMAGES} ${value}"
	    ;;
	--*)
	    usage "Unkown option ${key}"
	    ;;
	*)
	    ONLY_IMAGES="${ONLY_IMAGES} ${key}"
	    ;;
    esac
done

echo "Repository base is: ${DOCKER_REPO_BASE}"
echo "Build repository: ${BUILD_REPO}"
echo "Test repository: ${TEST_REPO}"

echo "Build images: ${BUILD_IMAGES}"
echo "Test images: ${TEST_IMAGES}"

if [ -n "${ONLY_IMAGES}" ]
then
    for include_image in ${ONLY_IMAGES}
    do
	image_exists=false
	for image in ${BUILD_IMAGES} ${TEST_IMAGES}
	do
	    if [ "${include_image}" = "${image}" ]
	    then
		image_exists=true
		break
	    fi
	done
	if [ "${image_exists}" = "false" ]
	then
	    usage "Unknown image specified: ${include_image}"
	fi
    done

    echo "Limiting to images: ${ONLY_IMAGES}"
fi

skip_image() {
    image=$1

    if [ -z "${ONLY_IMAGES}" ]
    then
	# All images, don't skip any.
	return 1
    fi
    for include_image in ${ONLY_IMAGES}
    do
	if [ "${image}" = "${include_image}" ]
	then
	    # Found it, don't skip this one.
	    return 1
	fi
    done
    # Not included, skip
    return 0
}

if [ -z "${DOCKER_REPO_BASE}" ]
then
    echo 'DOCKER_REPO_BASE must be set; perhaps run from bazel?' 1>&2
    exit 1
fi
BUILD_REPO="${DOCKER_REPO_BASE}/spiral/builder"
TEST_REPO="${DOCKER_REPO_BASE}/spiral/testbase"

cd tools/gendocker

# Using a private AWS repo? Login.
# aws ecr get-login-password --region us-west-2 |
#    docker login --username AWS --password-stdin ${DOCKER_REPO_BASE}

for IMAGE in ${BUILD_IMAGES}
do
    if skip_image ${IMAGE}
    then
	echo SKIPPING build image ${IMAGE}...
	continue
    fi
    echo BUILDING build image ${IMAGE}...
    tar -czh . |
	docker build -t ${BUILD_REPO}:${IMAGE} -t spiral/builder:${IMAGE} -f Dockerfile-${IMAGE}-build -
    if [ "${DO_PUSH}" = "true" ]
    then
	echo PUSHING ${IMAGE}
	docker push ${BUILD_REPO}:${IMAGE}
    fi
    echo DONE ${IMAGE}
done

for IMAGE in ${TEST_IMAGES}
do
    if skip_image ${IMAGE}
    then
	echo SKIPPING build image ${IMAGE}...
	continue
    fi
    echo BUILDING test image ${IMAGE}...
    tar -czh . |
	docker build -t ${TEST_REPO}:${IMAGE} -t spiral/testbase:${IMAGE} -f Dockerfile-${IMAGE} -
    if [ "${DO_PUSH}" = "true" ]
    then
	echo PUSHING ${IMAGE}
	docker push ${TEST_REPO}:${IMAGE}
    fi
    echo DONE ${IMAGE}
done

