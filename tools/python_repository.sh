#!/bin/bash

set -e

if [ $# -ne 1 ]
then
    echo "Usage: $0 <python version>" 1>& 2
    echo "Example: $0 3.6" 1>& 2
    exit 1
fi

VERSION=$1
SV=$(echo $1 | sed 's/\.//g')

for M in "m" ""
do
    DIR="/opt/python/cp${SV}-cp${SV}${M}"
    echo "Checking for ${DIR}"
    if [ -d ${DIR} ]
    then
	# Existing installation of python
	ln -sf ${DIR}/* .

	cat > BUILD <<EOF

# libpython target that doesn't link with libpython.  This should be used for
# external modules, since the python interpreter will already have libpython linked in.
cc_library(
  name="libpython",
  includes=["include/python${VERSION}${M}"],
  visibility=["//visibility:public"],
  hdrs=glob(["include/**/*.h"]),
)

py_runtime(
  name="runtime",
  files=[],
  interpreter="bin/python${VERSION}${M}",
  visibility=["//visibility:public"],
  python_version="PY3"
)

EOF
	exit 0
    fi
done

# Download and install python from APT
for PKG in python${VERSION}-minimal libpython${VERSION}-minimal libpython${VERSION}-dev
do
    apt-get download ${PKG}
    dpkg-deb -x ${PKG}*.deb .
done

cat > BUILD <<EOF

# libpython target that doesn't link with libpython.  This should be used for
# external modules, since the python interpreter will already have libpython linked in.
cc_library(
  name="libpython",
  includes=["usr/include/python${VERSION}", "usr/include"],
  visibility=["//visibility:public"],
  hdrs=glob(["usr/include/**/*.h"]),
)

py_runtime(
  name="runtime",
  files=[],
  interpreter="usr/bin/python${VERSION}",
  visibility=["//visibility:public"],
  python_version="PY3"
)

EOF

exit 0
