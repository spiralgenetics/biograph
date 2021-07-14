#!/usr/bin/env python3
""" Simple VCF exporter for BioGraph variants. """
#pylint: disable=mixed-indentation

from biograph import new_graph, reference, find_variants, version
from datetime import datetime
import sys
import argparse

PARSER = argparse.ArgumentParser(description="Export VCF entries for BioGraph variants",
    epilog="Example:\n\n$ vcfexport.py --in NA12878_S1.seqset --ref human_g1k_v37 --scaffold 2 --start 100000 --end 200000",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
PARSER.add_argument('--in', dest='filename', required=True, help='BioGraph file')
PARSER.add_argument('--ref', dest='reference', required=True, help='Reference directory')
PARSER.add_argument('--scaffold', '--chromosome', '--contig', dest='scaffold', required=True, type=str, help='Scaffold (chromosome) name')
PARSER.add_argument('--start', dest='start', required=True, type=int, help='Start position for assembly')
PARSER.add_argument('--end', dest='end', required=True, type=int, help='End position for assembly')

PARSER.add_argument('--min-size', dest='min_size', type=int, default=0, help='Minimum size of variant to report (default all)')
PARSER.add_argument('--min-overlap', dest='min_overlap', type=int, default=70, help='Minimum overlap (default 70)')
PARSER.add_argument('--version', action='version', version=version(), help='Show the libspiral version')
ARGS = PARSER.parse_args()

# Open the graph
bg = new_graph(ARGS.filename)
# Open the reference
ref = reference(ARGS.reference)
# Find variants in this range
all_variants = find_variants(bg, ref, ARGS.scaffold, ARGS.start, ARGS.end, ARGS.min_overlap)

if not all_variants:
    print "No variants found."
    sys.exit(1)

print """##fileformat=VCFv4.1
##fileDate={date}
##source=Spiral Genetics BioGraph v{version} {cmdline}
##reference=file://{reference}
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural Variant Type">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">
##INFO=<ID=AID,Number=.,Type=Integer,Description="Assembly IDs used in constructing this variant">
##INFO=<ID=AMBCOUNT,Number=1,Type=Integer,Description="Count of alternate locations for this end of an ambiguous breakend">
##INFO=<ID=AMBMATES,Number=1,Type=Integer,Description="Count of possible mate locations of an ambiguous breakend">
##INFO=<ID=SCL,Number=0,Type=Flag,Description="Possible SNP cluster">
##INFO=<ID=BEREF,Number=1,Type=Integer,Description="Average reference coverage across breakends">
##INFO=<ID=ENTROPYALT,Number=A,Type=Float,Description="Shannon entropy of alt allele if longer than 100 bp">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=TRANSPOSE,Number=1,Type=String,Description="Transposon FASTA sequence ID that this breakpoint anchor matches">
##INFO=<ID=SAS,Number=1,Type=Float,Description="Simple alignment score. Likelihood a breakend is not structural but rather aligns to reference simply">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Depth">
##FORMAT=<ID=ED,Number=1,Type=Integer,Description="Edit distance">
##FORMAT=<ID=OV,Number=A,Type=Integer,Description="Minimum read overlap in assembly">
##FORMAT=<ID=AQ,Number=A,Type=Integer,Description="Average quality">
##FILTER=<ID=low_depth,Description="The sample depth was below the minimum allowed threshold">
##FILTER=<ID=high_depth,Description="The sample depth was above the maximum allowed threshold">
##FILTER=<ID=low_edit_distance,Description="The edit distance is below the minimum allowed threshold">
##FILTER=<ID=too_many_alleles,Description="The set of possible alleles was too large/supported to be called">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">""".format(
    date=datetime.now().strftime('%Y%m%d'),
    version=version(),
    cmdline=" ".join(sys.argv[:]),
    reference=ARGS.reference
)

# TODO: needs genome::get_chromosome() for full scaffold size
#
# Since we only support one scaffold at a time, report length as the start of
# the first ref_range to the end of the last + offset.
#
print '##contig=<ID={id},length={length},url="file://{reffile}",offset={offset}>'.format(
    id=ARGS.scaffold,
    length=all_variants[-1].ref_range.end - all_variants[0].ref_range.start,
    reffile=ARGS.reference,
    offset=all_variants[0].ref_range.start
)

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"

for result in all_variants:
    # Filter variants by min-size
    filtered_variants = [v for v in result.variants if len(v.sequence) >= ARGS.min_size or v.range.size >= ARGS.min_size]

    if len(filtered_variants) == 0:
        # print "No variants passed filtering (out of {count} found).".format(count=len(variants))
        continue

    for i, v in enumerate(sorted(filtered_variants, key=lambda v: v.range.start), 1):
        # Must include the ref base prior to the event if no REF or ALT
        if not v.range.sequence or not v.sequence:
            position = v.range.start - 1
            last_base = str(ref.make_range(v.range.scaffold, position, position + 1, False).sequence)
            refseq = last_base + str(v.range.sequence)
            altseq = last_base + str(v.sequence)
        else:
            position = v.range.start
            refseq = v.range.sequence
            altseq = v.sequence

        # TODO: BLAST filtering
        filtered = "PASS"

        # TODO: group variants by position to call alleles for genotype
        if v.is_structural:
            gt = './.'
        else:
            gt = '0/1'

        # TODO: AD vs. DP.
        #
        # I think the reference part of AD should be:
        #
        #   rp = sum(result.coverage)/len(result.coverage)
        #
        # ...but this sometimes yields zero depth...?

        # Revisit when we implement allelic calling.

        # bool a0 = (depths[0] / tot_depth > .1);
        # bool a1 = (depths.size() > 1 ? depths[1] / tot_depth > .1 : false);
        # bool a2 = (depths.size() > 2 ? depths[2] / tot_depth > .1 : false);
        # if (call.alleles.size() > 3) { filter = "too_many_alleles"; }
        # else if (a0 && a1 && a2) { filter = "too_many_alleles";  }
        # else if (a0 && a1 && !a2) { gt = "0/1"; }
        # else if (a0 && !a1 && a2) { gt = "0/2"; }
        # else if (a0 && !a1 && !a2) { gt = "0/0"; }
        # else if (!a0 && a1 && a2) { gt = "1/2";  }
        # else if (!a0 && a1 && !a2) { gt = "1/1"; }
        # else if (!a0 && !a1 && a2) { gt = "2/2"; }

        # TODO: ED
        # Needs struct_var.cpp:sv_compute_edit_distance()

        if v.is_structural:
            print "{contig}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}\t{sample}".format(
                    contig=v.range.scaffold,
                    pos=position,
                    id=i,
                    ref=refseq if refseq else '.',
                    alt=altseq if altseq else '.',
                    qual='.',
                    filter=filtered,
                    info='NS=1;DP={dp};AID={aid}'.format(
                            dp=int(v.avg_depth),
                            aid=v.assembly_id
                        ),
                    fmt='GT:DP:OV',
                    sample='{gt}:{dp}:{ov}'.format(
                            gt=gt,
                            dp=int(v.avg_depth),
                            ov=v.min_overlap
                        )
                )
        else:
            print "{contig}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}\t{sample}".format(
                    contig=v.range.scaffold,
                    pos=position,
                    id=i,
                    ref=refseq if refseq else '.',
                    alt=altseq if altseq else '.',
                    qual='.',
                    filter=filtered,
                    info='NS=1;DP={dp};AID={aid}'.format(
                            dp=int(v.avg_depth),
                            aid=v.assembly_id,
                        ),
                    fmt='GT:DP:OV',
                    sample='{gt}:{dp}:{ov}'.format(
                            gt=gt,
                            dp=int(v.avg_depth),
                            ov=v.min_overlap
                        )
                )
