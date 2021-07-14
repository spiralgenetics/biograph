"""
VDB command wrapper
"""

import argparse
import sys

from biograph.vdb import (
    # vdb_delete,
    vcf_cmd,
    study_cmd,
    anno_cmd,
    query_cmd,
    # cohort,
    # trait,
    # phenotype
)

TOOLS = {
    "vcf": vcf_cmd.main,
    "study": study_cmd.main,
    "anno": anno_cmd.main,
    "query": query_cmd.main,
}

USAGE = """\
vdb - the Spiral Variant DataBase

    Subcommands:
        vcf              Import and export VCF files
        study            Gather, filter, and report on variants
        anno             Manage variant annotations

"""
    # Coming soon:
        # cohort           Group samples together into logical units
        # trait            Add or remove trait variables
        # phenotype        Assign trait values to samples
        # database         Manage VDB databases

    # secret:
        # query            Run arbitrary queries on the VDB

# database create, delete, list, show, copy

def parseArgs(clargs):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="vdb", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(clargs)

    TOOLS[args.cmd](args.options)

def main(clargs):
    ''' Access the variant database (beta) '''
    # NOTE: ^^^ this ^^^ docstring is used for the command description in __main__.py
    parseArgs(clargs)

if __name__ == '__main__':
    main(sys.argv[1:])
