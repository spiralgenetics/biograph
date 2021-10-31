# BioGraph

This is the codebase for the BioGraph genomic analysis platform. It is
released under the BSD 2-clause license. See the LICENSE file for details.

# Quick start

Download a pre-built binary from [GitHub releases](https://github.com/spiralgenetics/biograph/releases/).
Releases are statically built and should run on most recent Linux distributions.

In addition to the release, you will need a copy of the classifier model:

    wget https://archive.spiralgenetics.com/files/models/biograph_model-7.1.0.ml

or

    aws s3 cp s3://spiral-archive/models/biograph_model-7.1.0.ml . --request-payer requester

You will also need a few common bioinformatics packages. To install them on Ubuntu:

    sudo apt install -y vcftools tabix bcftools

BioGraph can also be run as a [Docker container](https://hub.docker.com/repository/docker/spiralgenetics/biograph):

    docker run spiralgenetics/biograph ...

# Documentation

[Full documentation on the GitHub wiki](https://github.com/spiralgenetics/biograph/wiki)

For other dependencies and full installation instructions, see the [online installation docs](https://github.com/spiralgenetics/biograph/wiki/Installation).

# Building the code

We recommend building on Ubuntu 18.04 or later. The following system dependencies are required:

    sudo apt-get install -y python3.6 python3.6-dev python3-distutils python3-apt python3-virtualenv virtualenv build-essential libbz2-dev libz-dev libcurl4-openssl-dev dh-autoreconf wamerican docker.io 

Note: Ubuntu 18.04.5 and later no longer include the python command by default. It is safest to avoid the system python and build from inside a virtualenv instead:

    virtualenv --python=python3.6 ~/buildenv
    . ~/buildenv/bin/activate
    
You will also need a recent version of [bazelisk](https://github.com/bazelbuild/bazelisk/releases/latest). 
Download the binary for your architecture and save it to a file named **bazel** somewhere in your PATH.

Build and run all tests:

    bazel test ... --config=manylinux2014

Create a release tarball:

    bazel build -c opt --config=release //python/biograph:package

The pip installable tarball will then be inside `bazel-bin/python/biograph/BioGraph*gz`
