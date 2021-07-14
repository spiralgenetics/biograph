#!/bin/bash
main() {
    set -e -x -o pipefail

    echo "fasta: '${fasta}'"
    echo "reference_name: '${reference_name}'"

    export REFDIR="reference/${reference_name}"

    echo "Downloading FASTA"
    mkdir -p fasta
    dx download "${fasta}" -o fasta/

    echo "Running biograph"
    mkdir -p reference
    biograph reference --in fasta/* --refdir ${REFDIR}

    echo "Archiving reference"
    tar zcf "${reference_name}.tgz" ${REFDIR}

    echo "Uploading reference archive"

    id=$(dx upload "${reference_name}.tgz" -p --path "reference/" --brief)
    dx-jobutil-add-output reference ${id} --class=file
}
