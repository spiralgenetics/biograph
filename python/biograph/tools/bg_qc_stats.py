"""
Basic QC Stats for a BioGraph file
"""

from __future__ import print_function

import sys
import random
import argparse
import math

import biograph


def map_reads(read_seq, mate_seq, ref):
    """
    Maps a pair of reads and returns their insert size
    Returns None if no exact match for the pairs
    """
    aln1 = ref.find(read_seq)

    if aln1.matches != 1:
        aln1 = ref.find(read_seq.rev_comp())
        if aln1.matches != 1:
            return None
        aln2 = ref.find(mate_seq)
    else:
        aln2 = ref.find(mate_seq.rev_comp())

    if aln2.matches != 1:
        return None

    read1 = aln1.get_match(0)
    read2 = aln2.get_match(0)

    if read1.chromosome != read2.chromosome:
        return None

    isize = max(read1.start, read1.end, read2.start, read2.end) - min(read1.start, read1.end, read2.start, read2.end)

    # if isize > 1000:
    #     print("huge insert size: {0} -> {1}, {2} -> {3} == {4}".format(read1.start, read1.end, read2.start, read2.end, isize))

    if isize < 0:
        raise RuntimeError(
            "negative insert size: {0} -> {1}, {2} -> {3} == {4}".format(read1.start, read1.end, read2.start, read2.end, isize))

    return isize


def estimate_insert(rm, ref, max_cnt=5000, sample_rate=0.90, outliers=0.05):
    """
    Estimate Mean, Median, SD of InsertSize of reads in a Readmap
    Outliers removes the top/bottom percent of insert sizes (typically chimeras, mis-mappings, or reference issues)
    """
    all_dists = []
    read_id = -1
    while read_id < rm.size() - 1 and max_cnt > 0:
        read_id += 1
        if random.random() >= sample_rate:
            continue

        read = rm.get_read_by_id(read_id)
        # Only want forward, paired reads - though the forward part might not matter
        if not read.has_mate():
            continue

        mate = read.get_mate()  # Runtime error, really?
        distance = map_reads(read.get_seqset_entry().sequence(),
                             mate.get_seqset_entry().sequence(),
                             ref)
        if distance is None:
            continue
        all_dists.append(distance)
        max_cnt -= 1

    all_dists.sort()
    obs = len(all_dists)

    # on unpaired data, obs == 0
    if obs == 0:
        return 0, 0, 0

    # trim off the top and bottom 5% edges outliers of lengths because they're
    # going to be from chimeric type seq
    trim = int(obs * outliers / 2)

    all_dists = all_dists[trim:-trim]

    obs = len(all_dists)
    mean = sum(all_dists) / float(obs)
    sd = math.sqrt(sum([math.pow(x - mean, 2) for x in all_dists]) / (obs - 1))
    return mean, all_dists[int(len(all_dists) / 2)], sd


def parse_args(args):
    """ Make pretty arguments """
    parser = argparse.ArgumentParser(prog="stats", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--biograph", metavar="BG", required=True,
                        help="BioGraph file containing an individual")
    parser.add_argument("-r", "--reference", metavar="REF",
                        help=("Reference genome folder. If not provided, insert size "
                              "and coverage are not estimated"))
    parser.add_argument("-s", "--sample", metavar="SAMPLE",
                        help="Accession id of sample to use. (default=all)")
    random.seed(1849282)
    args = parser.parse_args(args)
    return args


def print_stats(m_readmap, ref=None):
    """
    printing the stats
    """
    stats = m_readmap.get_pair_stats()
    print("{0: <19}{1:,}".format("NumReads:", stats.paired_reads + stats.unpaired_reads))
    bases = stats.paired_bases + stats.unpaired_bases
    print("{0: <19}{1:,}".format("NumBases:", bases))
    print("{0: <19}{1:,}".format("MaxReadLength:", m_readmap.max_read_len()))
    print("{0: <19}{1:,}".format("MinReadLength:", m_readmap.min_read_len()))
    print("{0: <19}{1:,}".format("NumPairedReads:", stats.paired_reads))
    print("{0: <19}{1:,}".format("NumUnpairedReads:", stats.unpaired_reads))
    print("{0: <19}{1:,}".format("NumPairedBases:", stats.paired_bases))
    print("{0: <19}{1:,}".format("NumUnpairedBases:", stats.unpaired_bases))

    if ref is not None:
        meanins, medianins, sdins = estimate_insert(m_readmap, ref)
        print("{0: <19}{1:.2f}".format("MeanInsertSize:", meanins))
        print("{0: <19}{1:.2f}".format("MedianInsertSize:", medianins))
        print("{0: <19}{1:.2f}".format("SDInsertSize:", sdins))
        print("{0: <19}{1:.2f}".format("EstimatedCoverage:", float(bases) / ref.size))
    print()

def main(args):
    ''' Get basic QC stats from a BioGraph '''
    # NOTE: ^^^ this ^^^ docstring is used for the command description in __main__.py
    args = parse_args(args)
    my_bg = biograph.BioGraph(args.biograph, biograph.CacheStrategy.MMAP)
    my_ref = None
    if args.reference is not None:
        my_ref = biograph.Reference(args.reference)
    if args.sample is not None:
        m_readmap = my_bg.open_readmap(args.sample)
        if m_readmap is not None:
            raise RuntimeError("Sample %s not found in BioGraph" % args.sample)
        print("{0: <19}{1:}".format("Sample:", args.sample))
        print_stats(m_readmap, my_ref)
    else:
        for sample in sorted(my_bg.metadata.samples):
            m_readmap = my_bg.open_readmap(sample)
            print("{0: <19}{1:}".format("Sample:", sample))
            print_stats(m_readmap, my_ref)

if __name__ == '__main__':
    main(sys.argv[1:])
