
import libspiral

FUZZ_DIST   = 1000   # Number of basepairs to fuzz
MIN_OVERLAP = 80     # Minimum overlap for assemblies
MAX_ANCHORS = 200    # Maximum anchors in fuzz region
MAX_STEPS   = 20000  # Maximum number of steps in graph walk

def refine_breakpoint(seqset, ref, contig, left, right):
    # Get reference ranges
    left_range = ref.make_range(contig, left - FUZZ_DIST, left + FUZZ_DIST, True)
    right_range = ref.make_range(contig, right - FUZZ_DIST, right + FUZZ_DIST, True)
    # Anchor reference
    left_anchors = libspiral.anchor(seqset, left_range, True, MIN_OVERLAP, MAX_ANCHORS)
    right_anchors = libspiral.anchor(seqset, right_range, False, MIN_OVERLAP, MAX_ANCHORS)
    # Perform assembly
    out = libspiral.assemble(left_anchors, right_anchors, MIN_OVERLAP, MAX_STEPS, True)
    # Filter out non-svs
    out_sv = [x for x in out if x.is_structural]
    # Error if more than one, return None is 0
    if len(out_sv) > 1:
        raise "Ambigious breakpoint"
    if len(out_sv) == 0:
        return None
    return out_sv[0]

def output_breakpoint(f, sv):
    print('\t'.join(map(str, [
        f, 
        sv.left_contig,
        sv.left_position,
        sv.right_position,
        sv.sequence,
        sv.min_overlap,
        sv.min_depth,
        sv.avg_depth,
        sv.assembly_sequence
    ])))

files = ["A00081_20150611_12_38", "A03226_2015034_14_10", "A02909_20150616_14_18", "A03614_2015034_14_10"]

locations = [
    ("4", 10211260, 10234567),
    ("13", 101894146, 101896525),
]

ref = libspiral.reference("data/baylor_hg19")
for f in files:
    seqset = libspiral.seqset("data/baylor_seqset/" + f + ".gbwt")
    for l in locations:
        sv = refine_breakpoint(seqset, ref, l[0], l[1], l[2])
        if sv:
            output_breakpoint(f, sv)


