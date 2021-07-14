"""
Program that will assemble over an arbitrary region of the genome

Picks anchors out of a reference region and attempts to reconstruct paths between
"""


from __future__ import print_function
import sys
import logging
import argparse
import functools

from heapq import heappush, heappop

import biograph

def parse_args(args):
    parser = argparse.ArgumentParser(prog="pair_tracer", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("biograph",
                        help="Biograph to use")
    parser.add_argument("--sample", default="",
                        help="Trace reads from the given sample")
    parser.add_argument("--min-insert-size", default=200, type=int,
                        help="Minimum insert size to search for (including the reads on both ends)")
    parser.add_argument("--max-insert-size", default=2000, type=int,
                        help="Maximum insert size to search for (including the reads on both ends)")
    parser.add_argument("--max-ambiguous-branches", default=2, type=int,
                        help="Prune search tree after encountering this many ambiguous branches")
    parser.add_argument("--min-overlap", default=30, type=int,
                        help="Do not trace through seqset entries smaller than this many bases")
    parser.add_argument("--debug", action="store_true",
                        help="Increase output verbosity")
    parser.add_argument("--max-steps", default=10000, type=int,
                        help="Maximum number of search steps to attempt before giving up")
    parser.add_argument("--direction", default=PAIRED_END, choices=[PAIRED_END, MATE_PAIR],
                        help=("Direction that the ends of pairs face "
                              "paired - end means the reads are oriented from the outside "
                              "in, and mate - pair means the reads are oriented from the "
                              "inside out"))
    args = parser.parse_args(args)
    setup_logging(args.debug)
    return args

def test_make_anchor():
    bg = biograph.BioGraph("/home/english/ajtrio/biographs/datasets/ajtrio/3.1.0/HG002.bg/")
    ref = biograph.Reference("/reference/hs37d5/")
    chrom, start, end = "1", 899922, 901922
    rm = bg.open_readmap()
    anchor_up = make_anchor(rm, ref, chrom, start)
    anchor_dn = make_anchor(rm, ref, chrom, end, False)
    eid = 0
    set_reads = set()
    rid = 0
    for i in anchor_up:
        print(">up_entry_%d\n%s" % (eid, i.entry.sequence()))
        for read in i.reads:
            set_reads.add(read[1])
            rid += 1
        eid += 1

    for rid, read in enumerate(set_reads):
        print(">up_read_%d\n%s" % (rid, read.get_seqset_entry().sequence()))

    eid = 0
    set_reads = set()
    rid = 0
    for i in anchor_dn:
        print(">dn_entry_%d\n%s" % (eid, i.entry.sequence()))
        for read in i.reads:
            set_reads.add(read[1])
            rid += 1
        eid += 1

    for rid, read in enumerate(set_reads):
        print(">dn_read_%d\n%s" % (rid, read.get_seqset_entry().sequence()))

class AnchorNode():
    """
    Queue-able and all that jazz
    """
    def __init__(self, entry, reads):
        self.entry = entry
        self.reads = reads

def make_anchor(m_rm, ref, chrom, position, upstream=True, anchor_size=70, max_window=300, max_anchors=5, m_anchors=None):
    """
    Given a readmap, reference, position, find a set of anchors' seqset entry
    if upstream, search 3'->5', else search downstream

    Works by returning all anchor_size kmers within max_window of the position that map uniquely within
    a reference and have a seqset entry - up to max_anchors are returned
    """
    if m_anchors is None:
        m_anchors = []

    if max_window == 0:
        return m_anchors

    #Make the anchor
    start, end = (position - anchor_size, position) if upstream else (position, position + anchor_size)
    try:
        anchor_seq = ref.make_range(chrom, start, end, True).sequence
    except RuntimeError as e:
        logging.warning("Unable to make_range on %s:%d - %s", chrom, position, str(e))
        return None

    position -= 10 if upstream else -50

    anchor_seq_rc = anchor_seq.rev_comp()

    # Ensure it's a unique anchor in reference
    ref_ctx_fwd = ref.find(anchor_seq)
    ref_ctx_rev = ref.find(anchor_seq_rc)
    matches = ref_ctx_fwd.matches + ref_ctx_rev.matches
    if matches != 1:
        return make_achor(m_rm, ref, chrom, position, upstream, anchor_size, max_window - 10, max_anchors, m_anchors)

    # Check to see if there's any coverage
    seq_entry = m_rm.seqset.find(anchor_seq)
    if not seq_entry:
        return make_achor(m_rm, ref, chrom, position, upstream, anchor_size, max_window -10, max_anchors, m_anchors)

    # Find reads starting here - or at least towards the beginning
    reads = [x for x in m_rm.get_reads_containing(seq_entry)]

    new_anchor = AnchorNode(seq_entry, reads)
    m_anchors.append(new_anchor)
    if len(m_anchors) >= max_anchors:
        return m_anchors

    #go find more
    return make_anchor(m_rm, ref, chrom, position, upstream, anchor_size, max_window - 1, max_anchors, m_anchors)


def tests():
    test_make_anchor()

if __name__ == '__main__':
    setup_logging(True)
    tests()
