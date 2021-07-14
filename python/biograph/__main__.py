#!/usr/bin/env python3
"""
biograph command wrapper
"""
import sys
import argparse

from biograph.tools import (
    bg_qc_stats,
    bgbinary_cmds,
    coverage,
    discover,
    export_aligned,
    install_tests,
    refhash,
    vdb,
)

from biograph.classifier import (
    qual_classifier,
    qual_classifier_PP
)

from biograph import version

# pylint: disable=unused-argument
def get_version(args):
    ''' Print the BioGraph version and exit '''
    print("biograph version %s" % version())

TOOLS = {
    "coverage": coverage.main,
    "create": bgbinary_cmds.create_cmd,
    "discovery":bgbinary_cmds.discovery_cmd,
    "exp_discover": discover.main,
    "export_aligned": export_aligned.main,
    "full_pipeline": bgbinary_cmds.full_pipeline,
    "install_test": install_tests.main,
    "qual_classifier": qual_classifier,
    "qual_classifier_PP": qual_classifier_PP,
    "reference": bgbinary_cmds.reference_cmd,
    "refhash": refhash.main,
    # "squareoff":bgbinary_cmds.squareoff,
    "stats": bg_qc_stats.main,
    "vdb": vdb.main,
    "version": get_version,
}

#squareoff       {TOOLS['squareoff'].__doc__}
USAGE = f"""\
biograph v{version()} - the BioGraph genome processing pipeline

    Pipeline Commands:
        full_pipeline   {TOOLS['full_pipeline'].__doc__}

        reference       {TOOLS['reference'].__doc__}
        create          {TOOLS['create'].__doc__}
        discovery       {TOOLS['discovery'].__doc__}
        coverage        {TOOLS['coverage'].__doc__}
        qual_classifier {TOOLS['qual_classifier'].__doc__}
        vdb             {TOOLS['vdb'].__doc__}

    Utility Commands:
        stats           {TOOLS['stats'].__doc__}
        version         {TOOLS['version'].__doc__}
        refhash         {TOOLS['refhash'].__doc__}

For help on any command, use the --help option:

    $ biograph full_pipeline --help

For full documentation: https://www.spiralgenetics.com/user-documentation

"""

ADDITIONAL = f"""\
    Additional commands:
        exp_discover    {TOOLS['exp_discover'].__doc__}
"""

def main():
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="biograph", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str,
                        help=argparse.SUPPRESS)
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help=argparse.SUPPRESS)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if '--help-all' in sys.argv:
        print(USAGE + ADDITIONAL, file=sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    TOOLS[args.cmd](args.options)
