"""
Looks at a project-level VCF and determines if paths of alleles create duplicate putative haplotypes.
The haplotype create from more VCF entries is removed. (beta)
"""
import os
import sys
import copy
import gzip
import tabix
import logging
import argparse
import functools

from StringIO import StringIO
from heapq import heappush, heappop
from collections import defaultdict, Counter

from biograph import BioGraph, Reference
from biograph.internal import GenomeGraph

@functools.total_ordering
class Path(object):

    """Holds a path through a GenomeGraph GNodes"""
    # pylint: disable=too-few-public-methods
    __slots__ = ['node', 'sequence', 'path']

    def __init__(self, node, sequence, path=None):
        # Node is the GenomeGraph Node that is at this point
        # sequence is the sequence we've built so far
        # path is a list of the previous GenomeGraph Nodes we've visited
        # seqset entry we can append a base to "sequence".
        self.node = node
        self.sequence = sequence
        self.path = path if path else [self.node]

    def extend_path(self, neigh, ref):
        """
        Return a new path that would be created from extending this path
        with the provided neighbor. ref is a BioGraph.Reference
        """
        new_seq = self.sequence + neigh.get_sequence(ref)
        new_path = list(self.path)
        new_path.append(neigh)
        return Path(neigh, new_seq, new_path)

    def __lt__(self, other):
        """
        Want to make sortable based self.node's downstream most position
        Need to fix the GNode's sorting, however
        """
        return self.node < other.node

    def __str__(self):
        return "<QueueNode len=%d seqlen=%d end=%s>" % (len(self.path), len(self.sequence), self.node)

    def __repr__(self):
        return str(self)


class VcfGraphBatcher(object):
    """
    Reads a vcf and creates a GenomeGraph in 'batches' within a buffer
    """

    def __init__(self, fn, ref, sub_chrom, sub_start, sub_end, distance=500):
        self.header = ""
        self.vcf_fn = fn
        self.ref = ref
        self.m_vcf = tabix.open(self.vcf_fn)
        self.all_data = set()
        self.sub_chrom = sub_chrom
        self.sub_start = sub_start
        self.sub_end = sub_end
        self.distance = distance

    def generate_graphs(self):
        """
        parses my region and yields each batch of
        chrom, min_start, max_end, GenomeGraph
        """
        def vcf_iter():
            """
            Yields batches of variants within distance
            """
            cur_batch = []
            prev_pos = None
            for data in self.m_vcf.query(self.sub_chrom, self.sub_start, self.sub_end):
                chrom = data[0]
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                end = pos + len(ref)
                if ',' in alt:
                    logging.error("VCF needs to have a single entry per line")
                    logging.error("Error line: %s" % "\t".join(data[:5]))
                    exit(1)

                anchor_trim = 1 if ref[0] == alt[0] else 0
                allele = alt[anchor_trim:]
                n_pos = pos + anchor_trim
                if cur_batch and prev_pos + self.distance < n_pos:
                    yield cur_batch
                    cur_batch = []

                prev_pos = max(prev_pos, n_pos)
                cur_batch.append((chrom, n_pos, end, allele, data))
            yield cur_batch

        for batch in vcf_iter():
            n_graph = GenomeGraph(self.ref)
            for var in batch:
                n_graph.add_var(*var)
            yield self.sub_chrom, batch[0][1], batch[-1][2], n_graph

def simple_dedup(m_graph, ref, t_chrom, t_start, t_end):
    """
    For each alt, make the haplotype sequence over the region
    If >1 alts make the same sequence, randomly pick one of them
    # This isn't exhaushative but whatever- pretty much truvari

    """
    #src = m_graph.get_node(t_chrom, t_start)
    #snk = m_graph.get_node(t_chrom, t_end + 1)
    src = m_graph.get_node(t_chrom, 0)
    snk = m_graph.get_node(t_chrom, ref.scaffold_lens[t_chrom])

    #start_seq = ref.make_range(src.chrom, src.abs_end - 1000, src.abs_end).sequence
    start_seq = ref.make_range(src.chrom, src.abs_end - 500, src.abs_end).sequence

    #end_seq = ref.make_range(src.chrom, max_end, max_end + 1000).sequence
    end_seq = ref.make_range(src.chrom, snk.abs_start, snk.abs_start + 500).sequence


    #src = m_graph.get_node(t_chrom, 0)
    #snk = m_graph.get_node(t_chrom, ref.scaffold_lens[t_chrom])
    # Have a problem with beginning/ends here.
    # like the t_end insnisn't coreect
    # Need to handle Reference regions out of range more gracefully
    #start_seq = ref.make_range(src.chrom, src.abs_end - 1000, src.abs_end).sequence
    #start_seq = ref.make_range(src.chrom, src.abs_end - 500, src.abs_end).sequence
    #end_seq = ref.make_range(src.chrom, max_end, max_end + 1000).sequence
    #end_seq = ref.make_range(src.chrom, snk.abs_start, snk.abs_start + 500).sequence

    def make_hap(alt_node):
        m_seq = StringIO()
        m_seq.write(str(start_seq))
        cur_node = src
        while cur_node != snk:
            if alt_node in m_graph.graph.successors(cur_node):
                m_seq.write(str(alt_node.get_sequence(ref)))
                cur_node = alt_node
                continue

            for edge in m_graph.graph.edges(cur_node):
                neighbor = edge[1]
                if not neighbor.is_alt and neighbor != alt_node:
                    m_seq.write(str(neighbor.get_sequence(ref)))
                    cur_node = neighbor

        m_seq.write(str(end_seq))
        m_seq.seek(0)
        return m_seq.read()

    logging.info("Creating putative haplotypes")
    all_haps = defaultdict(list)
    for alt_node in m_graph.alt_node_iter():
        if not alt_node.strand:
            continue
        m_hap = make_hap(alt_node)
        all_haps[m_hap].append(alt_node)

    for redund_seq in all_haps:
        redund = all_haps[redund_seq]
        if len(redund) > 1:
            for i in redund[1:]:
                m_graph.graph.remove_node(i)
    return m_graph

def smart_dedup(m_graph, ref, t_chrom, t_start, t_end):
    """
    Breadth first search.
    Keep a Queue of Nodes to trace.
    """
    src = m_graph.get_node(t_chrom, t_start)
    snk = m_graph.get_node(t_chrom, t_end + 1)

    # Need to handle Reference regions out of range more gracefully
    #start_seq = ref.make_range(src.chrom, src.abs_end - 1000, src.abs_end).sequence
    start_seq = ref.make_range(src.chrom, min(src.abs_end - 500, t_start), src.abs_end).sequence

    #end_seq = ref.make_range(src.chrom, max_end, max_end + 1000).sequence
    end_seq = ref.make_range(src.chrom, snk.abs_start, max(t_end, snk.abs_start + 500)).sequence

    logging.info("Creating putative haplotypes")

    path_heap = [Path(src, start_seq)]
    seen_seqs = defaultdict(list)
    seen_seqs[path_heap[0].sequence].append(path_heap[0])
    finished_vars = set()
    while path_heap:
        # get the first element off of the heap
        cur_path = heappop(path_heap)
        seen_seqs[cur_path.sequence].remove(cur_path)
        # For each edge
        for edge in m_graph.graph.edges(cur_path.node):
            neigh = edge[1]
            if neigh == snk:
                cur_path.node = neigh
                cur_path.path.append(neigh)
                cur_path.sequence = cur_path.sequence + end_seq
                new_path = cur_path
            else:
                new_path = cur_path.extend_path(neigh, ref)

            # ref node with shared sequence
            if not new_path.node.is_alt and len(seen_seqs[new_path.sequence]):
                ipos = 0
                while ipos < len(seen_seqs[new_path.sequence]):
                    other_path = seen_seqs[new_path.sequence][ipos]

                    if other_path.node != new_path.node:
                        ipos += 1
                        continue
                    #dedup elegible
                    seen_seqs[new_path.sequence].remove(other_path)
                    path_heap.remove(other_path)

                    o_alts = set([x for x in other_path.path if x.is_alt])
                    m_alts = set([x for x in new_path.path if x.is_alt])
                    a_set = o_alts - m_alts
                    b_set = m_alts - o_alts
                    if len(a_set) < len(b_set):
                        m_iter = b_set
                    else:
                        m_iter = a_set
                        # Gotta switch here for a_set because if
                        # this new_path is the 'non-parsimonous' one
                        # then I can't accidentally remove it's alts multiple time
                        new_path = other_path

                    for i in m_iter:
                        if i in m_graph.graph:
                            m_graph.graph.remove_node(i)

                        # Also what if it isn't in there

            if neigh != snk: # don't need to look aymore
                heappush(path_heap, new_path)
                seen_seqs[new_path.sequence].append(new_path)
    return m_graph

def create_haplotypes(m_graph, ref, t_chrom, t_start, t_end):
    src = m_graph.get_node(t_chrom, t_start)
    snk = m_graph.get_node(t_chrom, t_end + 1)

    #start_seq = ref.make_range(src.chrom, src.abs_end - 1000, src.abs_end).sequence
    start_seq = ref.make_range(src.chrom, min(src.abs_end - 500, t_start), src.abs_end).sequence

    #end_seq = ref.make_range(src.chrom, max_end, max_end + 1000).sequence
    end_seq = ref.make_range(src.chrom, snk.abs_start, max(t_end, snk.abs_start + 500)).sequence

    logging.info("Creating putative haplotypes")

    all_haps = defaultdict(list)
    redundant_haps = 0
    # Current slow part... change to StringIO

    for p in m_graph.all_paths(src, snk):
        cur_hap = StringIO(str(start_seq))
        for n in p[1:-1]:
            cur_hap.write(str(n.get_sequence(ref)))
        cur_hap.write(end_seq)
        cur_hap.seek(0)
        all_haps[cur_hap.read()].append(p)
        redundant_haps += 1

    logging.info("%d total haplotypes dedups to %d unique", redundant_haps, len(all_haps))
    return all_haps, redundant_haps != len(all_haps)


def dedup_entries(all_haps, all_vars):
    logging.info("Deduping VCF entries")
    dup_alts = set()
    for hap_path in all_haps.itervalues():
        if len(hap_path) == 1:  # unique - no variant to remove
            continue
        # I have more than one path to make a seq, let's look
        # at the difference between the two..
        # the first test they all had the same number of nodes.
        # and all these have 2 for redundant paths - so let's code for that
        a, b = set(hap_path[0]), set(hap_path[1])
        # add sets of redundant sets?
        dup_alts.add(frozenset([frozenset(a - b), frozenset(b - a)]))

    removed = 0
    for i in dup_alts:
        # Let's just get the one with the fewest number of vars
        max_alts = None
        for j in i:
            alts = [x for x in j if x.is_alt]
            if max_alts is None:
                max_alts = alts
            elif len(alts) > len(max_alts):
                max_alts = alts
        for j in max_alts:
            try:
                all_vars.remove(tuple(j.var_data))
                logging.debug("Removing %s", str(j.var_data))
                removed += 1
            except ValueError:
                # should log
                pass
    logging.info("Removed %d alleles", removed)
    return dup_alts


def main(args):
    args = parse_args(args)
    m_vcf = vcf_parser(args.variants)
    ref = Reference(args.reference)

    all_haps, any_redundancy = create_haplotypes(m_graph, ref, t_chrom, t_start, max_end)

def test_gbatcher():
    vcf = "/home/english/science/english/NovaSeq_HG02_replicates/FiveHG002.vcf.gz"
    ref = Reference("/reference/hs37d5/")
    chrom = "1"
    start = 669400
    end = 706368
    batcher = VcfGraphBatcher(vcf, ref, chrom, start, end)
    for chrom, start, end, graph in batcher.generate_graphs():
        print chrom, start, end
        for alt_node in graph.alt_node_iter():
            if alt_node.strand:
                print "\t".join(alt_node.var_data)
        print

def test_dedup():
    vcf = "/home/english/science/english/NovaSeq_HG02_replicates/FiveHG002.vcf.gz"
    vcf = "/home/english/science/english/NovaSeq_HG02_replicates/test22.vcf.gz"
    vcf = "/home/english/science/english/NovaSeq_HG02_replicates/test22.vcf.gz"
    vcf = "sid_chr22_input.vcf.gz"
    vcf = sys.argv[1]
    ref = Reference("/reference/hs37d5/")
    chrom = "3"
    start = 0
    end = ref.scaffold_lens[chrom]
    # start, end = 16242084, 16244015
    batcher = VcfGraphBatcher(vcf, ref, chrom, start, end)
    for chrom, start, end, graph in batcher.generate_graphs():
        logging.debug("Dedupping %s:%d-%d", chrom, start, end)
        try:
            #smart_dedup(graph, ref, chrom, start, end)
            simple_dedup(graph, ref, chrom, start, end)
        except RuntimeError:
            pass
        for alt_node in graph.alt_node_iter():
            if alt_node.strand:
                print "\t".join(alt_node.var_data)

        #for alt_node in graph.alt_node_iter():
            #if alt_node.strand:
                #all_vars.append(alt_node.var_data)

        #try:
            #all_haps, any_redundancy = create_haplotypes(graph, ref, chrom, start, end)
            #if any_redundancy:
                #duplicates = dedup_entries(all_haps, all_vars)
        #except RuntimeError:
            #pass

        #for entry in all_vars:
            #print "\t".join(entry)




if __name__ == '__main__':
    #test_gbatcher()
    setup_logging(True)
    test_dedup()

