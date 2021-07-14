#!/usr/bin/python3

"""
Demonstration program that shows how to trace between two seqset entries.

Picks random pairs out of a readmap and attempt to reconstruct the contents of the insert between them
"""

from __future__ import print_function
import sys
import random
import logging
import argparse
import functools

from heapq import heappush, heappop

import biograph

PAIRED_END = "paired-end"
MATE_PAIR = "mate-pair"


def parse_args(args):
    parser = argparse.ArgumentParser(prog="pair_tracer", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("biograph",
                        help="Biograph to use")
    parser.add_argument("--sample", default="",
                        help="Trace reads from the given sample")
    parser.add_argument("--count", default=10, type=int,
                        help="Number of pairs to attempt to trace")
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


def find_anchors(rm, ref, chrom, position, upstream=True, kmer_size=31, max_window=300):
    """
    give an position in the reference, find an anchor that's kmerSize that's
    unique within the reference. This iteratively works upstream or dnstream
    to find the anchor up max_window times. e.g. an upstream kmer will end upto
    max_window bp upstream from the provided position.
    Returns None if a unique hit was not found.
    """
    # need to check reference here
    if max_window < 1:  # no further to go
        return None
    if upstream:
        try:
            anchor_ref = ref.make_range(chrom, position - kmer_size, position, True)
        except RuntimeError as e:
            logging.warning("Unable to make_range on %s:%d - %s", chrom, position, str(e))
            return None
        position -= 1
    else:
        try:
            anchor_ref = ref.make_range(chrom, position, position + kmer_size, True)
        except RuntimeError as e:
            logging.warning("Unable to make_range on %s:%d - %s", chrom, position, str(e))
            return None
        position += 1
    ref_ctx = ref.find(anchor_ref.sequence)
    if ref_ctx.matches != 1 or len(my_bg.read_search(anchor_ref.sequence)) == 0:
        return make_anchor(my_bg, ref, chrom, position, upstream, kmer_size, max_window - 1)
    return anchor_ref


def trace_anchors(m_rm, args):
    """
    Close to the old anchored assembly approach.
    Given a reference and start/end coordinates, find the first
    location where a while read maps exactly upstream/downstream.
    Use those as the two pairs and send to the same method as
    trace anchors
    """
    up_anchor = find_anchor()
    dn_anchor = find_anchor()

    tracer = Tracer(rm, args)
    result = tracer.trace_between_entries(read, mate)
    if result:
        logging.info("Found path:")
        logging.info("%s", result)
    else:
        logging.info("No path found")
    return result


def trace_random_pair(m_rm, args):
    """
    Picks a random pair out of the readmap and attempts to trace between the two ends
    """
    # Find a read that has a mate
    while True:
        read_id = random.randint(0, m_rm.size() - 1)
        read = m_rm.get_read_by_id(read_id)
        if read.has_mate():
            break
    # If it's not in is original orientation, get the original orientation
    if not read.is_original_orientation():
        read = read.get_rev_comp()
    assert read.is_original_orientation()

    # Find the other end of the pair.
    mate = read.get_mate()

    # If we're using "mate-pair" instead of "paired-end", we have to
    # trace the opposite direction to find the mate.
    if args.direction == MATE_PAIR:
        read = read.get_rev_comp()
        mate = mate.get_rev_comp()

    logging.info("Tracing between:")
    logging.info("Pair: %s", read.get_seqset_entry().sequence())
    logging.info("Mate: %s",  mate.get_seqset_entry().sequence())

    tracer = Tracer(rm, args)
    result = tracer.trace_between_reads(read, mate)
    if result:
        logging.info("Found path:")
        logging.info("%s", result)
    else:
        logging.info("No path found")
    return result


class Tracer(object):

    """Tracer holds the state for an individual trace between two seqset entries"""
    DIR_READ = True
    DIR_MATE = False
    DIR_NAMES = {DIR_READ: "READ", DIR_MATE: "MATE"}
    MAX_ASM_LENGTH = 1000
    MIN_ASM_LENGTH = 200
    MAX_STEPS = 10000
    MIN_OVERLAP = 70
    MAX_AMBIG = 2

    @functools.total_ordering
    class QueueEntry(object):

        """QueueEntry represents a path in the search queue"""
        # pylint: disable=too-few-public-methods
        __slots__ = ['entry', 'direction', 'sequence', 'ambiguous_count']

        def __init__(self, entry, direction, sequence, ambiguous_count):
            # Seqset entry that this entry is searching from.  This
            # seqset entry is facing the opposite direction of
            # "sequence", so when we use "push_front_drop" on the
            # seqset entry we can append a base to "sequence".
            self.entry = entry
            self.direction = direction
            self.sequence = sequence
            self.ambiguous_count = ambiguous_count

        def __lt__(self, other):
            self_ambig = self.ambiguous_count
            other_ambig = other.ambiguous_count
            if self_ambig != other_ambig:
                return self_ambig < other_ambig
            return len(self.sequence) < len(other.sequence)

        def __str__(self):
            return "<QueueEntry len=%d entrylen=%d direction=%s seq=%s ambig=%d>" % (
                len(self.sequence),
                len(self.entry),
                Tracer.DIR_NAMES[self.direction],
                self.sequence,
                self.ambiguous_count)

    def __init__(self, m_readmap):
        self.m_readmap = m_readmap
        self.trace_ends = {Tracer.DIR_READ: {}, Tracer.DIR_MATE: {}}
        self.queue = []
        self.seen = set()

    def add_to_queue(self, entry, direction, sequence, ambiguous_count):
        """
        Checks to see if an entry matches a trace from the other direction.

        If the given entry matches a trace from the other direction, returns the total sequence
        of the trace.  Otherwise adds a new entry to the search queue and returns None.
        """
        for read in self.m_readmap.get_prefix_reads(entry):
            rc_mate_id = read.get_rev_comp().get_read_id()
            other_end_sequence = self.trace_ends[not direction].get(rc_mate_id, None)
            if other_end_sequence:
                assert sequence[-len(read):] == other_end_sequence[-len(read):].rev_comp()
                total_sequence = sequence[
                    :-len(read)] + other_end_sequence.rev_comp()
                logging.debug("Found a candidate path of length %d:", len(total_sequence))
                logging.debug("This end:  %s", sequence)
                logging.debug("Other end: %*s%s", len(sequence) - len(read), "", other_end_sequence.rev_comp())
                if (len(total_sequence) <= Tracer.MAX_ASM_LENGTH and
                        len(total_sequence) >= Tracer.MIN_ASM_LENGTH):
                    # Found one!
                    return total_sequence
                else:
                    logging.debug("Candidate path exceeded insert size limits")
            self.trace_ends[direction][read.get_read_id()] = sequence
        heappush(self.queue, self.QueueEntry(entry, direction, sequence, ambiguous_count))

    def trace_between_entries(self, up_entry, dn_entry):
        """
        Attempts to find a path between a read and its mate.

        Returns the total sequence of the path if successful.
        """
        self.add_to_queue(up_entry, Tracer.DIR_READ, up_entry.sequence(), 0)
        self.add_to_queue(dn_entry, Tracer.DIR_MATE, dn_entry.sequence(), 0)
        assert len(self.queue) == 2

        step_count = 0
        while self.queue:
            step_count += 1
            if step_count > Tracer.MAX_STEPS:
                logging.debug("Step count %d exceeded maximum with %d items left in the queue",
                              step_count, len(self.queue))
                return
            queue_entry = heappop(self.queue)
            result = self.trace_next(queue_entry)
            if result is not None:
                # Successfully found a path.
                logging.debug("Path found after %d steps", step_count)
                return result

    def trace_between_reads(self, read, mate):
        """
        Attempts to find a path between a read and its mate.

        Returns the total sequence of the path if successful.
        """
        # Should only be called once per Tracer object
        assert not self.queue

        read_seq = read.get_seqset_entry().sequence()
        mate_seq = mate.get_seqset_entry().sequence()

        # "read" and "mate" have their ends pointing towards each
        # other.  Our primary operation is "push_front", so when
        # traversing the seqset we need to trace through using the
        # reverse complement.
        read_entry = read.get_rev_comp().get_seqset_entry()
        mate_entry = mate.get_rev_comp().get_seqset_entry()

        return self.trace_between_entries(read_entry, mate_entry)

    def trace_next(self, qentry):
        """
        Given an entry in the search queue, trace all possible paths forward from it.
        """
        logging.debug("Processing search queue entry: %s", qentry)
        if len(qentry.sequence) > Tracer.MAX_ASM_LENGTH:
            logging.debug("Sequence length exceeded maximum insert size")
            return
        next_entries = {}
        next_lens = []
        for base, comp_base in zip("ACGT", "TGCA"):
            entry = qentry.entry.push_front_drop(comp_base, Tracer.MIN_OVERLAP)
            if not entry:
                logging.debug("Adding next base %s violates minimum overlap; skipping", base)
                continue
            drop_count = len(qentry.entry) + 1 - len(entry)
            if drop_count:
                logging.debug("Adding next base %s drops %d bases from seqset entry", base, drop_count)
            next_entries[base] = entry
            next_lens.append(len(entry))
        next_lens.sort()

        if len(next_lens) > 1:
            # Any path forward that has <= the length of the second
            # best path forward we mark as ambiguous and don't search
            # unless we can't find a way forward otherwise.
            ambiguous_len = next_lens[-2]
        else:
            ambiguous_len = 0

        for next_base, next_entry in next_entries.items():
            next_ambiguous_count = qentry.ambiguous_count
            if len(next_entry) <= ambiguous_len:
                # Length of this entry is <= the second best entry, so
                # count it as ambiguous
                next_ambiguous_count += 1
                if next_ambiguous_count > Tracer.MAX_AMBIG:
                    logging.debug("Adding base %s exceeds ambiguous count", next_base)
                    continue
                logging.debug("Adding base %s with an ambiguous path", next_base)
            # else:
            #    logging.debug("Adding base %s", next_base)
            new_seq = qentry.sequence + next_base
            logging.debug("Full sequence: %s", new_seq)
            result = self.add_to_queue(next_entry, qentry.direction, new_seq, next_ambiguous_count)
            if result is not None:
                return result

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    bg = biograph.BioGraph(args.biograph)
    rm = bg.open_readmap(args.sample)

    tot_found, tot_not_found = 0, 0
    for i in range(args.count):
        if trace_random_pair(rm, args):
            tot_found += 1
        else:
            tot_not_found += 1

    logging.info("%d total traces successful; %d unsuccessful", tot_found, tot_not_found)
