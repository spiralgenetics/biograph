# BioGraph

This is the codebase for the BioGraph genomic analysis platform. It is free for academic use. See the LICENSE file for details.

# Quick start

Download a pre-built binary from [GitHub releases](https://github.com/spiralgenetics/biograph/releases/).
Releases are statically built and should run on most recent Linux distributions.

In addition to the release, you will need a copy of the classifier model:

    $ wget https://spiral-public.s3.us-west-2.amazonaws.com/models/biograph_model-7.0.0.ml

or

    $ aws s3 cp s3://spiral-public/models/biograph_model-7.0.0.ml . --no-sign-request

You will also need a few common bioinformatics packages. To install them on Ubuntu:

    $ sudo apt install -y vcftools tabix bcftools

BioGraph can also be run as a [Docker container](https://hub.docker.com/repository/docker/spiralgenetics/biograph):

    $ docker run spiralgenetics/biograph ...

# Documentation

For other dependencies and full installation instructions, see the [online installation docs](https://www.notion.so/spiralgenetics/Installation-8105bf74808a4b9e822a77fd4ad6cbf2).

# Building the code

We recommend building on Ubuntu 18.04 or later. The following system dependencies are required:

    $ apt-get install -y python3.6 python3.6-dev build-essential libbz2-dev libz-dev libcurl4-openssl-dev dh-autoreconf wamerican docker.io

You will also need a recent version of [bazelisk](https://github.com/bazelbuild/bazelisk/releases/latest). 
Download the binary for your architecture and save it to a file named **bazel** somewhere in your PATH.

Build and run all tests:

    $ bazel test ... --config=manylinux2014

Create a release tarball:

    $ bazel build -c opt --config=release //python/biograph:package
