#!/usr/bin/env python3
"""
validate.py

Validate calls in a VCF file against a sequence set of the same individual.
Adds a BD INFO field to each VCF entry indicating the minumum reference and
allelic coverage.

Only variants with a filter value of PASS are considered. Structural variants
are not yet supported.
"""
from __future__ import print_function
import sys
import argparse

import vcf

from biograph import Seqset, Readmap, Reference, Sequence

# pylint: disable=too-many-locals


def validate(seq, read_map, ref, infile, outfile, margin=100):
    ''' Run the validation '''

    # load the seqset and readmap
    print("Loading Seqset...")
    my_sample = Seqset(seq)
    print("Loading Readmap...")
    my_readmap = Readmap(read_map)

    # this reference should match the reference used in the input VCF
    print("Loading reference...")
    ref = Reference(ref)

    # input and output VCF files
    input_vcf = open(infile, 'r')
    output_vcf = open(outfile, 'w')

    # set up the VCF reader
    vcf_reader = vcf.Reader(input_vcf)

    # add a new INFO field for our output VCF
    # pylint: disable=protected-access
    vcf_reader.infos.update(
        {
            'BD': vcf.parser._Format(
                id='BD', num=None, type='Integer',
                desc='Minimum BioGraph depths found for reference and alt alleles'
            )
        }
    )

    # set up the VCF writer
    # Note: this writes out the header immediately, so make any necessary header changes first.
    vcf_writer = vcf.Writer(output_vcf, vcf_reader)

    validated = 0
    skipped = 0

    print("Processing", infile)
    # walk the VCF
    for variant in vcf_reader:
        # Pass through variants with a non-empty FILTER, SVs, or insufficent context
        if variant.FILTER or variant.var_type == 'sv' or variant.POS < margin:
            skipped += 1
            vcf_writer.write_record(variant)
            continue

        try:
            # sequence before this position
            pre = str(
                ref.make_range(
                    variant.CHROM,
                    variant.POS - margin,
                    variant.POS,
                    False
                ).sequence
            )

            # sequence after the reference at this position
            post = str(
                ref.make_range(
                    variant.CHROM,
                    variant.POS + len(variant.REF),
                    variant.POS + len(variant.REF) + margin, False
                ).sequence
            )

        except RuntimeError:
            # Too close to the end of the contig, so pass it through.
            skipped += 1
            vcf_writer.write_record(variant)
            continue

        # reference sequence for the whole region
        seq = pre + variant.REF + post
        # find the reference coverage for this region
        ref_cov = my_sample.seq_coverage(seq, my_readmap)

        # start one base before the variant
        first = margin - 1
        # end one base after the affected reference region
        last = margin + len(variant.REF) + 1
        # first BD entry is the minimum reference coverage
        bdent = [min(ref_cov[first:last])]

        # compute BD for each allele
        for allele in variant.ALT:
            # build the alternate sequence for the whole region
            newseq = pre + str(allele.sequence) + post
            # end one base after the allele
            last = margin + len(allele.sequence) + 1
            # compute the coverage for this allele
            seq_cov = my_sample.seq_coverage(newseq, my_readmap)
            # compute the coverage for the reverse complement
            rev_cov = my_sample.seq_coverage(Sequence(newseq).rev_comp(), my_readmap)[::-1]
            # update BD with the minimum allele coverage
            bdent.append(min([x + y for x, y in zip(seq_cov[first:last], rev_cov[first:last])]))

        # Write out the VCF record for this variant
        variant.add_info('BD', bdent)
        vcf_writer.write_record(variant)
        validated += 1

        if not validated % 10000:
            print(validated, "records validated")

    output_vcf.close()
    print('validated:', validated, 'skipped:', skipped)


def run(args):
    """
    Executing script
    """
    parser = argparse.ArgumentParser(
        description="BioGraph VCF validator\n\n"
        "Validate calls in a VCF file against a sequence set of the same individual.\n"
        "Adds a BD INFO field to each VCF entry indicating the minumum reference\n"
        "and allelic coverage.\n"
        "\n"
        "Only variants with a filter value of PASS are considered. Structural variants\n"
        "are not yet supported.\n",
        epilog="Example:\n\n$ validate.py --seqset NA12878.seqset --readmap NA12878.readmap --ref human_g1k_v37 --in my_calls.vcf --out validated.vcf\n",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--seqset', dest='seqset', required=True,
                        help='Sequence set file of the individual used for the input VCF')
    parser.add_argument('--readmap', dest='readmap', required=True,
                        help='Readmap of the individual used for the input VCF')
    parser.add_argument('--ref', dest='reference', required=True,
                        help='Reference directory. Must match reference used for input VCF')
    parser.add_argument('--in', dest='infile', required=True, help='Input VCF file')
    parser.add_argument('--out', dest='outfile', required=True, help='Output VCF file')
    parser.add_argument('--margin', dest='margin', type=int, default=100,
                        help='Number of bases to inspect on either side of each variant (default: 100)')
    my_args = parser.parse_args(args)

    validate(my_args.seqset, my_args.readmap, my_args.reference, my_args.infile, my_args.outfile, my_args.margin)

# Run as a simple command line tool
if __name__ == '__main__':
    run(sys.argv[1:])
