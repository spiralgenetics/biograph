#!/bin/sh
# Install biograph and run install tests.  Called from inside a container.

echo INSTALL TEST running for image ${RULENAME}

set -e

. /home/spiral/venv/bin/activate

CACHEOPT=""
if [ -d /mnt/.cache ]
then
    mkdir -p /mnt/.cache/shared-pip-http
    mkdir -p /mnt/.cache/pip-${RULENAME}
    if ! test -e /mnt/.cache/pip-${RULENAME}/http
    then
	ln -s /mnt/.cache/shared-pip-http /mnt/.cache/pip-${RULENAME}/http
    fi
    CACHEOPT="--cache-dir=/mnt/.cache/pip-${RULENAME}"
fi

# Try to avoid pulling src packages if possible
pip install --prefer-binary ${CACHEOPT} /dist/${PIP_PACKAGE_FN}

# Make sure the model is cached in memory since otherwise reading it
# could cause the install test to time out.
echo Caching ${MODEL} from /share...
cat /share/releases/${MODEL} > /dev/null

/home/spiral/venv/bin/biograph install_test -m /share/releases/${MODEL}
