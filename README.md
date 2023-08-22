# BioGraph

This is the codebase for the BioGraph genomic analysis platform. It is
released under the BSD 2-clause license. See the LICENSE file for details.

# Quick start

Download a pre-built binary from [GitHub releases](https://github.com/spiralgenetics/biograph/releases/).
Releases are statically built and should run on most recent Linux distributions.

In addition to the release, you will need a copy of the classifier model available on [Zenodo](https://zenodo.org/record/8273311)

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

# Additional Resources

[BioGraph Primer](https://medium.com/@robflickenger/here-be-dragons-gaps-in-the-human-genome-and-how-inverting-the-analysis-can-help-e6a8078168ed) A primer on the challenges in population scale bioinformatics and the design choices behind BioGraph.

[White Paper Publication 2020](https://www.biorxiv.org/content/10.1101/2020.04.24.060202v1) Leveraging a WGS compression and indexing format with dynamic graph references to call structural variants.

[ASHG 2019](https://archive.spiralgenetics.com/files/www/ASHG2019Poster.pdf) A graph genome of the Arab population with base pair accurate structural variation for population genotyping. A collaboration with SIDRA Medical.

[AGBT 2017](https://archive.spiralgenetics.com/files/www/AGBT-2017-Poster.pdf) BioGraph: a reference-agnostic and rapidly queryable NGS read data format enabling flexible analysis at scale. A collaboration with University of Texas Health Science Center and Baylor College of Medicine

[HUGO 2016](https://archive.spiralgenetics.com/files/www/ReducingMendelianinconsistenciesinatriodatasetwitharapidqueryformat.pdf) Using read overlap assembly to accurately identify structural genetic differences in an Ashkenazi Jewish trio. A collaboration with the Baylor College of Medicine Human Genome Sequencing Center.

[AGBT 2016](https://archive.spiralgenetics.com/files/www/BioGraphPosterAGBT2016.pdf) BioGraph Suite: A format to analyze multiple graph genomes created directly from read data

