#!/usr/bin/env python3

"""
Demonstration program that shows how to trace between two seqset entries.

trace_through_pairs picks random pairs from a readmap and attempts to trace through the seqset to reconstruct the insert between them.

"""

from __future__ import print_function
import argparse
from heapq import heappush, heappop
import functools
import random
import typing

import biograph

PAIRED_END = "paired-end"
MATE_PAIR = "mate-pair"

parser = argparse.ArgumentParser(
    description="Pick random pairs out of a readmap and attempt to reconstruct the contents of the insert between them")
parser.add_argument("biograph", help="Biograph to use")
parser.add_argument(
    "--readmap", default="", help="Open the given readmap")
parser.add_argument(
    "--count", default=10, help="Number of pairs to attempt to trace")
parser.add_argument(
    "--min-insert-size", default=200, help="Minimum insert size to search for (including the reads on both ends)")
parser.add_argument(
    "--max-insert-size", default=2000, help="Maximum insert size to search for (including the reads on both ends)")
parser.add_argument(
    "--max-ambiguous-branches", default=2, help="Prune search tree after encountering this many ambiguous branches")
parser.add_argument(
    "--min-overlap", default=30, help="Do not trace through seqset entries smaller than this many bases")
parser.add_argument(
    "-v", "--verbose", action="count", help="Increase output verbosity")
parser.add_argument(
    "--max-steps", default=10000, help="Maximum number of search steps to attempt before giving up")
parser.add_argument(
    "--direction", default=PAIRED_END, choices=[PAIRED_END, MATE_PAIR],
    help="Direction that the ends of pairs face; paired-end means the reads are oriented from the outside in, and mate-pair means the reads are oriented from the inside out")

args = parser.parse_args()

bg = biograph.BioGraph(args.biograph)
rm = bg.open_readmap(args.readmap)

def trace_random_pair():
    """Picks a random pair out of the readmap and attempts to trace between the two ends"""
    # Find a read that has a mate
    while True:
        read_id = random.randint(0, rm.size() - 1)
        read = rm.get_read_by_id(read_id)
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

    print("Tracing between:")
    print("Pair: %s" % read.get_seqset_entry().sequence())
    print("Mate: %s" % mate.get_seqset_entry().sequence())

    tracer = Tracer()
    result = tracer.trace_between(read, mate)
    if result:
        print("Found path:")
        print("%s" % result)
    else:
        print("No path found")
    print("")
    return result


class Tracer(object):
    """Tracer holds the state for an individual trace between two seqset entries"""
    DIR_READ = True
    DIR_MATE = False
    DIR_NAMES = {DIR_READ: "READ", DIR_MATE: "MATE"}

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

    def __init__(self):
        self.trace_ends = {Tracer.DIR_READ: {}, Tracer.DIR_MATE: {}}
        self.queue = []
        self.seen = set()

    def add_to_queue(self, entry, direction, sequence, ambiguous_count):
        """Checks to see if an entry matches a trace from the other direction.

If the given entry matches a trace from the other direction, returns the total sequence
of the trace.  Otherwise adds a new entry to the search queue and returns None.
"""
        for read in rm.get_prefix_reads(entry):
            rc_mate_id = read.get_rev_comp().get_read_id()
            other_end_sequence = self.trace_ends[
                not direction].get(rc_mate_id, None)
            if other_end_sequence:
                assert sequence[-len(read):] == other_end_sequence[
                    -len(read):].rev_comp()
                total_sequence = sequence[
                    :-len(read)] + other_end_sequence.rev_comp()
                if args.verbose:
                    print(
                        "Found a candidate path of length %d:" %
                        len(total_sequence))
                    print("This end:  %s" % sequence)
                    print(
                        "Other end: %*s%s" %
                        (len(sequence) - len(read), "", other_end_sequence.rev_comp()))
                if (len(total_sequence) <= args.max_insert_size and
                        len(total_sequence) >= args.min_insert_size):
                    # Found one!
                    return total_sequence
                else:
                    if args.verbose:
                        print("Candidate path exceeded insert size limits")
            self.trace_ends[direction][read.get_read_id()] = sequence
        heappush(self.queue, self.QueueEntry(
            entry, direction, sequence, ambiguous_count))

    def trace_between(self, read, mate):
        """Attempts to find a path between a read and its mate.

Returns the total sequence of the path if successful."""
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

        self.add_to_queue(
            read_entry, Tracer.DIR_READ, read_seq, 0)
        self.add_to_queue(
            mate_entry, Tracer.DIR_MATE, mate_seq, 0)
        assert len(self.queue) == 2

        step_count = 0
        while self.queue:
            step_count += 1
            if step_count > args.max_steps:
                if args.verbose >= 1:
                    print(
                        "Step count %d exceeded maximum with %d items left in the queue" %
                        (step_count, len(self.queue)))
                return
            queue_entry = heappop(self.queue)
            result = self.trace_next(queue_entry)
            if result is not None:
                # Successfully found a path.
                if args.verbose:
                    print("Path found after %d steps" % step_count)
                return result

    def trace_next(self, qentry):
        """Given an entry in the search queue, trace all possible paths forward from it."""
        if args.verbose >= 2:
            print("Processing search queue entry: %s" % qentry)
        if len(qentry.sequence) > args.max_insert_size:
            if args.verbose >= 2:
                print("Sequence length exceeded maximum insert size")
            return
        next_entries = {}
        next_lens = []
        for base, comp_base in zip("ACGT", "TGCA"):
            entry = qentry.entry.push_front_drop(comp_base, args.min_overlap)
            if not entry:
                if args.verbose >= 3:
                    print(
                        "Adding next base %s violates minimum overlap; skipping")
                continue
            if args.verbose >= 3:
                drop_count = len(qentry.entry) + 1 - len(entry)
                if drop_count:
                    print(
                        "Adding next base %s drops %d bases from seqset entry" %
                        (base, drop_count))
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
                if next_ambiguous_count > args.max_ambiguous_branches:
                    if args.verbose >= 2:
                        print(
                            "Adding base %s exceeds ambiguous count" %
                            next_base)
                    continue
                if args.verbose >= 2:
                    print("Adding base %s with an ambiguous path" % next_base)
            else:
                if args.verbose >= 2:
                    print("Adding base %s" % next_base)
            new_seq = qentry.sequence + next_base
            if args.verbose >= 3:
                print("Full sequence: %s" % new_seq)
            result = self.add_to_queue(
                next_entry, qentry.direction, new_seq, next_ambiguous_count)
            if result is not None:
                return result

tot_found, tot_not_found = 0, 0
for i in range(args.count):
    if trace_random_pair():
        tot_found += 1
    else:
        tot_not_found += 1

print("%d total traces successful; %d unsuccessful" % (tot_found, tot_not_found))
