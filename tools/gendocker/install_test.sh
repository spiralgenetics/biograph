#!/bin/bash

set -e

RULENAME="$1"
DOCKER_IMAGE="$2"
LAMBDA_TOY_DATA="$3"
PIP_PACKAGE="$4"
MODEL="$5"

PIP_PACKAGE_FN=$(basename ${PIP_PACKAGE})

cat > Dockerfile-installtest-${RULENAME} <<EOF
FROM ${DOCKER_IMAGE}

USER root
ADD --chown=spiral:spiral ${LAMBDA_TOY_DATA} /

RUN if [ ! -d /lambdaToyData ]; then ln -s /lambdaToyData* /lambdaToyData; fi
WORKDIR /lambdaToyData
RUN mkdir /dist
COPY ${PIP_PACKAGE_FN} /dist/${PIP_PACKAGE_FN}
COPY run_install_test.sh /dist/run_install_test.sh

USER spiral
RUN mkdir -p /home/spiral/.cache
EOF

tar -czh Dockerfile-installtest-${RULENAME} \
    -C ${PWD}/tools/gendocker run_install_test.sh \
    -C /share/releases ${LAMBDA_TOY_DATA} \
    -C ${PWD}/$(dirname ${PIP_PACKAGE}) ${PIP_PACKAGE_FN} | \
    docker build \
	   -f Dockerfile-installtest-${RULENAME} \
	   -t installtest-${RULENAME} \
	   -

CACHEBIND=""
if [ -d /mnt/.cache ]
then
    mkdir -p /mnt/.cache/shared-pip-http
    chmod 1777 /mnt/.cache/shared-pip-http
    mkdir -p /mnt/.cache/pip-${RULENAME}
    chmod 1777 /mnt/.cache/pip-${RULENAME}
    CACHEBIND="--mount type=bind,destination=/mnt/.cache,src=/mnt/.cache"
fi

docker run \
       --mount type=bind,destination=/share,src=/share,readonly \
       ${CACHEBIND} \
       --env "RULENAME=${RULENAME}" \
       --env "PIP_PACKAGE_FN=${PIP_PACKAGE_FN}" \
       --env "TEST_CMD=${TEST_CMD}" \
       --env "MODEL=${MODEL}" \
       installtest-${RULENAME}:latest \
       /dist/run_install_test.sh
