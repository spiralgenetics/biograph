#!/bin/bash
main() {
    set -e -x -o pipefail

    echo "input_reads: '${input_reads}'"
    echo "reference: '${reference}'"
    echo "sample: '${sample_id}'"
    echo "create_options: '${create_options}'"
    echo "variants_options: '${variants_options}'"

    export BIOGRAPH="${sample_id}.bg"

    echo "Downloading reference"
    dx download "$reference" -o - | tar zxf -

    export REFDIR=$(echo reference/*)

    echo "Running biograph"
    dx download "$input_reads" -o - | biograph create --in - --out ${BIOGRAPH} --ref ${REFDIR} --tmp . ${create_options}

    echo "Calling variants"
    biograph variants --in ${BIOGRAPH} --out "${BIOGRAPH}/variants.vcf" --ref ${REFDIR} ${variants_options}

    vcf-sort < "${BIOGRAPH}/variants.vcf" | bgzip > "${BIOGRAPH}/variants.vcf.gz"
    tabix -p vcf "${BIOGRAPH}/variants.vcf.gz"

    echo "Uploading BioGraph"
    for f in metadata/bg_info.json qc/create_log.txt qc/create_stats.json qc/kmer_quality_report.html qc/variants_log.txt qc/variants_stats.json seqset; do
        id=$(dx upload ${BIOGRAPH}/$f -p --path "${BIOGRAPH}/$f" --brief)
        dx-jobutil-add-output $(basename ${f/.*}) ${id} --class=file
    done

    id=$(dx upload ${BIOGRAPH}/variants.vcf.gz -p --path "${BIOGRAPH}/" --brief)
    dx-jobutil-add-output variants ${id} --class=file

    id=$(dx upload ${BIOGRAPH}/variants.vcf.gz.tbi -p --path "${BIOGRAPH}/" --brief)
    dx-jobutil-add-output variants_index ${id} --class=file

    id=$(dx upload ${BIOGRAPH}/coverage/*.readmap -p --path "${BIOGRAPH}/coverage/" --brief)
    dx-jobutil-add-output readmap ${id} --class=file
}
