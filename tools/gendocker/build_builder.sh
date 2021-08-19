#!/bin/sh

# Generate docker images needed to build biograph

set -e
set -v

# cd tools/gendocker

tar -czh . |
    docker build -t spiralgenetics/manylinux2014:noref -f Dockerfile-manylinux2014-build-noref -


tar -czh . |
    docker build -t spiralgenetics/manylinux2014:latest -f Dockerfile-manylinux2014-build -

