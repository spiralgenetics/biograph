#!/usr/bin/python3

import argparse
from libspiral import *

parser = argparse.ArgumentParser(description='Genotype some folks.')
parser.add_argument('--ref', help='Which refernce to use')
parser.add_argument('--dels', help='The deletion list as a .bed file')
parser.add_argument('seqsets', metavar='seqset', nargs='+',
                   help='SEQSETs for each individual')

args = parser.parse_args()

print "Opening reference"
ref = reference(args.ref)
with open(args.dels, 'r') as content_file:
    content = content_file.read()
dels = [x.split('\t') for x in content.split('\n')][:-1]

print "Loading SEQSETs"
seqsets = {}
read_len = 0
for f in args.seqsets:
    seqsets[f] = seqset(f)
    if read_len == 0:
        read_len = seqsets[f].read_len
    if read_len != seqsets[f].read_len:
        raise "Mismatched read lengths"

min_overlap = int(.7*read_len)

def seek_dels(gb, d):
    print "Seeking", d
    slop = int(3 * read_len)
    max_anchors = slop / 2
    max_steps = 10 * slop
    readset = gb.all_readset()
    left = int(d[1])
    right = int(d[2])
    middle = int((left + right) / 2)
    left_range = ref.make_range(d[0], left - slop, left + slop, False)
    right_range = ref.make_range(d[0], right - slop, right + slop, False)
    left_anchor = find_anchors(gb, left_range, True, min_overlap, max_anchors, readset)
    right_anchor = find_anchors(gb, right_range, False, min_overlap, max_anchors, readset)
    results = assemble(left_anchor, right_anchor, min_overlap, max_steps, True, readset)
    if len(results) == 0:
        return []
    return [x for x in results[0] if x.is_structural]

del_vars = []
for fn in args.seqsets:
    print "Computing variation in", fn
    del_vars += [[seek_dels(seqsets[fn], d) for d in dels]]
    
print "Zippering and flattening"
zip_vars = zip(*del_vars)
flat_vars = [[x for y in z for x in y] for z in zip_vars]

def normalize_asm(var):
    if var.assembly_begin < min_overlap:
        return ""
    if var.assembly_end + min_overlap > len(var.assembly_sequence):
        return ""
    if not var.left_forward or not var.right_forward:
        return ""
    if var.left_position >= var.right_position:
        return ""
    return str(var.assembly_sequence)[var.assembly_begin - min_overlap : var.assembly_end + min_overlap]

def normalize_list(l):
    if len(l) == 0:
        return []
    if len(l) > 2:
        return []  # This is questionable
    alts = {}
    ref_contig = l[0].left_contig
    ref_start = None
    ref_end = 'X'
    for v in l:
        vf = v.flip()
        asm = normalize_asm(vf)
        if asm == "":
            continue
        alts[asm] = vf
        ref_start = max(ref_start, vf.left_position)
        ref_end = min(ref_end, vf.right_position)
    if ref_start == None:
        return []
    rr = ref.make_range(ref_contig, ref_start, ref_end, True)
    alleles = [rr]
    for k, v in alts.iteritems():
        alleles += [[sequence(k), v]]
    return alleles

print "Normalizing"
var_sets = [normalize_list(x) for x in flat_vars]

def genotype(gb, alls):
    if len(alls) == 0:
        return [0,0]
    ref_cov = gb.coverage(alls[0], gb.all_readset())
    alleles = [sum(ref_cov)/len(ref_cov)]
    for x in alls[1:]:
        alt_cov = gb.seq_coverage(x[0], gb.all_readset())
        alleles += [min(alt_cov[50:-50])]
    bins = [int(bool(x)) for x in alleles]
    if sum(bins) == 0 or sum(bins) > 2:
        return None
    both = []
    for i in range(len(bins)):
        if bins[i]:
            both += [i]
    if len(both) == 1:
        both = [both[0], both[0]]
    return both

print "Genotyping"

for i in range(len(var_sets)):
    vs = var_sets[i]
    d = dels[i]
    line = ""
    if len(vs) == 0:
        line += "%s\t%d\t%d\t" % (d[0], int(d[1]), int(d[2]))
    else:
        v = vs[1][1]
        line += "%s\t%d\t%d\t" % (v.left_contig, v.left_position, v.right_position-1)
        #bits = []
        #for a in vs[1:]:
        #    bits += [""]
        #line += " ".join(bits) + "\t"
    for fn in args.seqsets:
        gt = genotype(seqsets[fn], vs)
        if gt is None:
            line += "?/?\t"
        else:
            line += "%d:%d\t" % (gt[0], gt[1])
    print line
    
    
