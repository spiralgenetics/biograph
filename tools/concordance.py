#!/usr/bin/env python3
''' GIAB concordance for the human weekly run. '''

# pylint: disable=missing-docstring

from __future__ import print_function

import re
import subprocess

# s3cmd --skip-existing get s3://spiral-repo/giab/SNPs.vcf.gz GIAB_SNPs.vcf.gz
# s3cmd --skip-existing get s3://spiral-repo/giab/indels.vcf.gz GIAB_indels.vcf.gz

# s3cmd --skip-existing get s3://spiral-repo/giab/variants.bed.gz GIAB_variants.bed.gz
# gunzip -c GIAB_variants.bed.gz > GIAB_variants.bed

def snp_concordance():
    # download GIAB SNPs
    subprocess.call(['s3cmd', '--skip-existing', 'get',
        's3://spiral-repo/giab/SNPs.vcf.gz',
        'GIAB_SNPs.vcf.gz'
    ])

    # isolate anchored assembly SNPs
    subprocess.call(['vcftools',
        '--vcf', 'nonSV_bed_filtered_primitives.vcf',
        '--remove-indels',
        '--recode', '--recode-INFO-all',
        '--out', 'SNPs'
    ])

    # compare anchored assembly SNPs with GIAB SNPs
    subprocess.call(['vcftools',
        '--vcf', 'SNPs.recode.vcf',
        '--gzdiff', 'GIAB_SNPs.vcf.gz',
        '--out', 'SNP_concordance'
    ])

def indel_concordance():
    # download GIAB indels
    subprocess.call(['s3cmd', '--skip-existing', 'get',
        's3://spiral-repo/giab/indels.vcf.gz',
        'GIAB_indels.vcf.gz'
    ])

    # isolate anchored assembly indels
    subprocess.call(['vcftools',
        '--vcf', 'nonSV_bed_filtered_primitives.vcf',
        '--keep-only-indels',
        '--recode', '--recode-INFO-all',
        '--out', 'indels'
    ])

    # compare anchored assembly indels with GIAB indels
    subprocess.call(['vcftools',
        '--vcf', 'indels.recode.vcf',
        '--gzdiff', 'GIAB_indels.vcf.gz',
        '--out', 'indel_concordance'
    ])

def prepare():
    subprocess.call(['s3cmd', '--skip-existing', 'get',
        's3://spiral-repo/giab/variants.bed.gz',
        'GIAB_variants.bed.gz'
    ])

    subprocess.call(['s3cmd', '--skip-existing', 'get',
        's3://spiral-repo/tools/info_filter.py',
    ])

    subprocess.call(['s3cmd', '--skip-existing', 'get',
        's3://spiral-repo/tools/vcfallelicprimitives',
    ])

    subprocess.call(['chmod', '+x', 'vcfallelicprimitives'])

    with open('GIAB_variants.bed', 'w') as fout:
        subprocess.call(['gunzip', '-c', 'GIAB_variants.bed.gz'], stdout=fout)

    # Sort the VCF output from Anchored
    with open('variants.sorted.vcf', 'w') as fout:
        subprocess.call(['vcf-sort', 'variants.vcf'], stdout=fout)

    # isolate the non-SVs
    with open('nonSV.vcf', 'w') as fout:
        subprocess.call(['vcf_filter.py',
            '--no-filtered',
            '--local-script', 'info_filter',
            'variants.sorted.vcf', 'sv'
        ], stdout=fout)

    # Filter the anchored assembly output vcf file of nonSVs with the GIAB bed file
    subprocess.call(['vcftools',
        '--vcf', 'nonSV.vcf',
        '--bed', 'GIAB_variants.bed',
        '--recode', '--recode-INFO-all',
        '--out', 'nonSV_bed_filtered'
    ])

    # break up alt alleles representing multiple mismatches
    # into the simplest possible separate records
    with open('nonSV_bed_filtered_primitives.vcf', 'w') as fout:
        subprocess.call(['./vcfallelicprimitives', 'nonSV_bed_filtered.recode.vcf'], stdout=fout)

    snp_concordance()
    indel_concordance()

FORMAT_HEAD = '{:10} {:<15} {:<15} {:<15} {:<15} {:<15}'
FORMAT_ITEM = '{:10} {:<15} {:<15} {:<15} {:<15%} {:<15%}'

def analyze(name, filename):
    pattern = re.compile('Found ([0-9]+) SNPs (.*) file')
    with open(filename, 'r') as fin:
        for line in fin:
            match = pattern.search(line)
            if match:
                value = match.group(1)
                part = match.group(2)
                if 'common' in part:
                    common = float(value)
                elif 'main' in part:
                    anchored_only = float(value)
                elif 'second' in part:
                    giab_only = float(value)

    sensitivity = common / (common + giab_only)
    ppv = 1 - (anchored_only / (common + anchored_only))
    print(FORMAT_ITEM.format(name, int(common), int(anchored_only), int(giab_only), sensitivity, ppv))

def main():
    prepare()

    print(FORMAT_HEAD.format('NAME', 'COMMON', 'ANCHORED', 'GIAB', 'SENSITIVITY', 'PPV'))
    analyze('SNP', 'SNP_concordance.log')
    analyze('indel', 'indel_concordance.log')


if __name__ == "__main__":
    main()
