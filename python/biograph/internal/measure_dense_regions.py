"""
Given an input VCF

Count for every SV (>=50bp), the number of neighboring SVs within BUFF

Optionally writes to CNT file each SVs:
    chrom   start_pos   num_neighbor_svs

For those CNTs, we then create bed regions spanning calls with more than DEN neighbors
and write that bed to OUT

The OUT.bed can then be `bedtools merge`'d  to reduce redundancy in the regions.
Then the VCF can have all entries within OUT.bed removed with `bedtools subtract`

The OUT.bed's 4th column is the average_number_of_neighbors,number_of_events
"""
import sys
import gzip
import bisect
import argparse

from collections import Counter, defaultdict

def parse_args(clargs):
    parser = argparse.ArgumentParser(prog="measure_dense_regions", description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", metavar="VCF", required=True,
                        help="Input to process")
    parser.add_argument("-b", "--buffer", metavar="BUFF", type=int, default=1000,
                        help="Span of the regions to count neighbors (%(default)s)")
    parser.add_argument("-d", "--density", metavar="DEN", type=int, default=100,
                        help="Create bed of regions with at least (%(default)s) neighbors")
    parser.add_argument("-c", "--counts", metavar="CNT", default=None,
                        help="Write each SV's number of neighbors to file")
    parser.add_argument("-o", "--output", metavar="OUT", default=None,
                        help="Output bed file name (stdout)")
    args = parser.parse_args(clargs)
    return args

def get_type_lens(entry):
    """
    Parse an entry and return it's sv_type and it's sv_len
    """
    mREF = entry[3]
    # TODO - should get the longest?
    mALTs = entry[4].split(',')
    sv_types = []
    sv_lens = []
    # Get type for counting - MYTYPES
    for mALT in mALTs:
        if len(mREF) == len(mALT):
            sv_types.append("REPL")
            sv_lens.append(len(mREF))
        elif len(mREF) == 1:
            sv_types.append("INS")
            sv_lens.append(len(mALT) - 1)
        elif len(mALT) == 1:
            sv_types.append("DEL")
            sv_lens.append(len(mREF) - 1)
        elif len(mREF) > len(mALT):
            sv_types.append("SUBSDEL")
            sv_lens.append(len(mREF) - len(mALT))
        elif len(mALT) > len(mREF):
            sv_types.append("SUBSINS")
            sv_lens.append(len(mALT) - len(mREF))
        else:
            logging.error(str(entry))
            logging.error("shouldn't have some new crazy type\n")
            exit()

    # MYSIZES
    ret_lens = []
    for sv_len in sv_lens:
        if sv_len < 10:
            ret_lens.append("1-9")
        elif sv_len < 50:
            ret_lens.append("10-49")
        elif sv_len < 300:
            ret_lens.append("50-299")
        elif sv_len < 1000:
            ret_lens.append("300-999")
        else:
            ret_lens.append("gt1000")
    
    return sv_types, ret_lens

def main(args):
    args = parse_args(args)
    
    fh = gzip.GzipFile(args.input, 'r') if args.input.endswith(".gz") else open(args.input, 'r')
    bed_out = open(args.output, 'w') if args.output is not None else sys.stdout
    if args.counts is not None:
        cnt_out = open(args.counts, 'w')
    
    svpositions = defaultdict(list)
    for line in fh:
        if line.startswith("#"):
            continue
        data = line.strip().split('\t')
        ty, ln = get_type_lens(data)
        if len([x for x in ln if x not in ["1-9", "10-49"]]):
            svpositions[data[0]].append(int(data[1]))
    
    cnt = Counter()
    for chrom in svpositions:
        cur_bed_entry = [None, None, None]
        sum_density = 0
        cnt_density = 0.0
        data = sorted(svpositions[chrom])
        for cur_pos in data:
            start = bisect.bisect_left(data, cur_pos - args.buffer)
            end = bisect.bisect_right(data, cur_pos + args.buffer)
            num_neighbors = end - start - 1
            
            if args.counts is not None:
                cnt_out.write("%s\t%d\t%d\n" % (chrom, cur_pos, num_neighbors))
            
            if num_neighbors >= args.density:
                if cur_bed_entry[0] is None:
                    cur_bed_entry = [chrom, cur_pos, cur_pos + 1]
                else:
                    cur_bed_entry[2] = cur_pos
                sum_density += num_neighbors
                cnt_density += 1
            else:
                if cur_bed_entry[0] is not None:
                    bed_out.write("%s\t%d\t%d\t%.2f,%d\n" % tuple(cur_bed_entry + [sum_density/cnt_density, cnt_density]))
                cur_bed_entry = [None, None, None]
                sum_density = 0.0
                cnt_density = 0
        
    if cur_bed_entry[0] is not None: # close it up
        bed_out.write("%s\t%d\t%d\t%.2f,%d\n" % tuple(cur_bed_entry + [sum_density/cnt_density, cnt_density]))

if __name__ == '__main__':
    main(sys.argv[1:])
