#!/usr/bin/env python3
# encoding=utf8
"""
    visualize.py

    Print a text visualization of a given assembly.
"""
from __future__ import print_function
import sys
import argparse


def coverage_print(var, region, use_ascii=True):
    """
    Prints coverage only
    """
    my_seq = var.ref_range.sequence
    prev = -1
    did_dots = False
    p_vert = chr(0x2502)
    for pos, cov in enumerate(var.coverage[:-1]):
        if cov != prev:
            print(u"{:<12}:{:>12}   {:2} {:.1}{}".format(region.chrom.decode(), region.start + pos, min(99, cov),
                                                         my_seq[pos], "|" if use_ascii else p_vert).encode('utf-8'))
            prev = cov
            did_dots = False
        elif not did_dots:
            sys.stdout.write((" " * 32) + '.\n')
            did_dots = True
    print(u"{:<12}:{:>12}   {:2} {:.1}{}".format(region.chrom, region.start + len(var.coverage),
                                                 min(99, var.coverage[-1]), my_seq[-1], "|"
                                                 if use_ascii else p_vert).encode('utf-8'))


def main(args):
    """
    Executing script
    """
    # pylint: disable=too-many-branches
    import biograph as bgsdk
    from biograph.utils import load_regions
    parser = argparse.ArgumentParser(
        description="Print a text visualization for a given region of a Seqset.\n"
        "\n"
        "The Readmap parameter is optional, but improves coverage reporting when available.\n",
        epilog="Example:\n\n$ biograph visualize --biograph NA12878.bg --ref human_g1k_v37/ --region 2:195000-198000",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--biograph", required=True, help="BioGraph file to visualize")
    parser.add_argument("-r", "--ref", dest="reference", required=True, help="Reference directory")
    parser.add_argument("-R", "--region_bed", type=str, help="Scaffold (a.k.a. chromosome or contig) name")
    parser.add_argument(
        "--region", type=str, nargs="*", help="Region over which to analyze variants in form chrom:start-end.")
    parser.add_argument("--min-size", type=int, default=0, help="Minimum size of variant to report (default all)")
    parser.add_argument("--min-overlap", dest="min_overlap", type=int, default=70, help="Minimum overlap (default 70)")
    parser.add_argument("--ascii", dest="use_ascii", action="store_true", help="Use ASCII instead of UTF8")
    parser.add_argument(
        "--breakpoint", action="store_true", help="Use find_breakpoint_variants instead of find_region_variants")
    parser.add_argument(
        "--buffer",
        default=300,
        type=int,
        help="Buffer argument for find_breakpoint_variants or amount to ± from region")
    parser.add_argument("-c", "--coverage", action="store_true", help="Print coverage if no variants found")
    parser.add_argument(
        "-s", "--sample", default=None, type=str, help="Name of sample when running with merged BioGraph")
    parser.add_argument("--version", action="version", version=bgsdk.version(), help="Show the BioGraph version")
    parser.add_argument("--debug", action="store_true", help="Show additional debug information on error")
    my_args = parser.parse_args(args)

    if not my_args.debug:
        sys.tracebacklimit = 0

    # Open the graph
    try:
        my_bg = bgsdk.BioGraph(my_args.biograph)
    except RuntimeError:
        raise RuntimeError('Could not open BioGraph ' + my_args.biograph)
    if len(my_bg.metadata.samples) != 1:
        if my_args.sample is None:
            raise RuntimeError("Merged BioGraph needs --sample to be specified")
        if my_args.sample not in my_bg.metadata.samples:
            raise RuntimeError("Sample %s not found in BioGraph. Available: %s" % (my_args.sample, ", ".join(
                my_bg.metadata.samples.keys())))
        my_bg.set_readmap(my_args.sample)

    # Open the reference
    try:
        ref = bgsdk.Reference(my_args.reference)
    except RuntimeError:
        raise RuntimeError('Could not open reference ' + my_args.reference)

    # Find variants for this reference range
    for region in load_regions(my_args.region_bed, my_args.region):
        if region.chrom not in ref.chromosomes:
            raise RuntimeError('Requested chromosome ' + region.chrom + ' does not exist in reference.')

        if my_args.breakpoint:
            all_variants = bgsdk.find_breakpoint_variants(
                my_bg,
                ref=ref,
                supercontig=region.chrom,
                start=region.start,
                end=region.end,
                min_overlap=my_args.min_overlap,
                buf_len=my_args.buffer)
        else:
            all_variants = bgsdk.find_region_variants(
                my_bg,
                ref=ref,
                supercontig=region.chrom,
                start=max(0, region.start - my_args.buffer),
                end=region.end + my_args.buffer,
                min_overlap=my_args.min_overlap)

        sys.stdout.write(str(region) + '\n')

        if not isinstance(all_variants, list):
            all_variants = [all_variants]

        for var in all_variants:
            if (not var.variants) and my_args.coverage:
                coverage_print(var, region, my_args.use_ascii)
            else:
                bgsdk.visualize(var, min_size=my_args.min_size, use_ascii=my_args.use_ascii)


if __name__ == "__main__":
    main(sys.argv[1:])
