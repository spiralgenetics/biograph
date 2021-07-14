"""
Given an BioGraph discovery vcf, reduce the graph complexity by removing extra 'noise' in SNP/INDELS that
have PDP==0 and aren't in phase with SVs. Inputs and outputs MUST be sorted.
"""
import sys
import argparse

from collections import defaultdict

import pysam

def parse_args(args):
    """ Make pretty arguments """
    parser = argparse.ArgumentParser(prog="rvg_cleaner", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--variants", metavar="VCF", default="/dev/stdin",
                        help="Input VCF file to parse (%(default)s)")
    parser.add_argument("-o", "--output", metavar="OUT", default="/dev/stdout",
                        help="Output VCF to write (%(default)s)")
    args = parser.parse_args(args)
    return args


def main(args):
    """
    Main
    """
    args = parse_args(args)
    # We can split this by chromosome UNTIL we get inter-chromosomal translocations going
    m_vcf = pysam.VariantFile(args.variants)
    output = pysam.VariantFile(args.output, 'w', header=m_vcf.header.copy()) # pylint: disable=no-member
    backlog = defaultdict(list) # PI: [list of variants we're waiting on]
    known_inphase_lookup = {} # PIs of things to keep
    cur_chrom = None
    for entry in m_vcf:
        if cur_chrom != entry.chrom:
            del(backlog)
            del(known_inphase_lookup)
            backlog = defaultdict(list)
            known_inphase_lookup = {}
            cur_chrom = entry.chrom
        m_pi = entry.samples[0]["PI"]
        if "SVLEN" in entry.info:
            output.write(entry)
            # put PI in list of 'always outputs'
            known_inphase_lookup[m_pi] = True
            # flush backlog
            for bent in backlog[m_pi]:
                output.write(bent)
            try:
                # not guaranteed to have a backlog
                del(backlog[m_pi])
            except KeyError:
                pass
        elif entry.samples[0]["PDP"] != 0 or m_pi in known_inphase_lookup:
            output.write(entry)
        else:
            backlog[m_pi].append(entry)

if __name__ == '__main__':
    main(sys.argv[1:])
