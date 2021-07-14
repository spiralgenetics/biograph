"""
Runs all the steps needed to create a project-level VCF
"""
import sys
import json
import logging
import argparse

import biograph
import biograph.tools as bgtools
from biograph.utils import bgtool_args

def parse_args(args):
    """"
    #parse args
    """
    parser = argparse.ArgumentParser(prog="pvcf_pipe", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #only really need the pcmp parameters
    #so add a group around that
    #multiple biographs - positional at the end
    parser.add_argument(
        "-b", "--biograph", metavar="BG", required=True, help="Merged BioGraph file containing individuals")
    parser.add_argument("-r", "--reference", metavar="REF",
                        required=True, help="Reference genome folder")
    parser.add_argument("-o", "--output", metavar="OUT",
                        default="output.vcf", help="Annotated vcf file output.")
    parser.add_argument("-p", "--pedigree", default=None,
                        help="Pedigree file for samples in the biographs")

    #pedigree file is required - only samples in the pedigree file
    # will be processed
    #Really the only parameters we need, currently
    parser.add_argument(
        "-m", "--min-insert", default=200, type=int, help="Minimum insert size to consider paired (%(default)s)")
    parser.add_argument(
        "-M", "--max-insert", default=1000, type=int, help="Maximum insert size to consider paired (%(default)s)")
    parser.add_argument(
        "-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true", help="Verbose logging")
    # and a temp dir
    args = parser.parse_args(args)
    setup_logging()
    return args

def main(args):
    args = parse_args(args)
    # lowered bar
    """
    For each biograph - sample_lookup = {bg'ssample: path to this bg}
    ensure every sample is accounted for in the pedigree file.
       otherwise exit
    For each biograph - vcfs = [vcf_file in the biograph]
    make - temp_merge file
    run merge_vcfs vcfs > temp_merge_file
    For sample in sample_lookup:
        pcmp_lookup = {sample:temp_pcmp_output}
        run pcmp on temp_merge_file > pcmp_lookup[temp_pcmp_output]
    run vcf_sample_paste pcmp_lookup.values() | meanno > args.output
    """
    # just a series of tempfile.mkstemp and bgtools.thing.run("formatted from what I got")

if __name__ == '__main__':
    main(sys.argv[1:])

