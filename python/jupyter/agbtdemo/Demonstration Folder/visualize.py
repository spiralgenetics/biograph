#!/usr/bin/env python3
""" Text-based visualizer for BioGraph variants. """
from biograph import new_graph, reference, find_variants, version, visualize
import sys
import argparse

PARSER = argparse.ArgumentParser(description="Visualize BioGraph variants",
	epilog="Example:\n\n$ visualize.py --in NA12878_S1.seqset --ref human_g1k_v37 --scaffold 2 --start 100000 --end 200000",
	formatter_class=argparse.RawDescriptionHelpFormatter
)
PARSER.add_argument('--in', dest='filename', required=True, help='BioGraph file to visualize')
PARSER.add_argument('--ref', dest='reference', required=True, help='Reference directory')
PARSER.add_argument('--scaffold', '--chromosome', '--contig', dest='scaffold', required=True, type=str, help='Scaffold (chromosome) name')
PARSER.add_argument('--start', dest='start', required=True, type=int, help='Start position for assembly')
PARSER.add_argument('--end', dest='end', required=True, type=int, help='End position for assembly')
# PARSER.add_argument('--readset', dest='readset', default=None, help='Readset bitmap to use (default: all readsets)')
PARSER.add_argument('--min-size', dest='min_size', type=int, default=0, help='Minimum size of variant to report (default all)')
PARSER.add_argument('--min-overlap', dest='min_overlap', type=int, default=70, help='Minimum overlap (default 70)')
PARSER.add_argument('--ascii', action='store_true', default=False, help='Use ASCII instead of UTF8')
PARSER.add_argument('--version', action='version', version=version(), help='Show the libspiral version')
ARGS = PARSER.parse_args()

# Open the graph
bg = new_graph(ARGS.filename)
# Open the reference
ref = reference(ARGS.reference)
# Make a reference range for the area of interest
ref_range = ref.make_range(ARGS.scaffold, ARGS.start, ARGS.end, False)

# if ARGS.readset:
# 	readset = bg.load_readset(ARGS.readset)
# else:
# 	readset = bg.all_readset()

# Find variants for this reference range
all_variants = find_variants(bg, ref, ARGS.scaffold, ARGS.start, ARGS.end, ARGS.min_overlap)

if len(all_variants) == 0:
	print "No variants found."
	sys.exit(1)

for v in all_variants:
	visualize(v.variants, v.coverage, v.ref_range, min_size=ARGS.min_size, ascii=ARGS.ascii)
