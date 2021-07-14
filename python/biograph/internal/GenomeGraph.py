"""
Represents a Genome as a directed graph capable of adding variants as new nodes and edges
"""

from __future__ import print_function

import bisect
import logging
import random
from itertools import islice
from io import BytesIO
from collections import defaultdict

import networkx as nx
from biograph import Sequence

import vcf
STRREVCMP = str.maketrans("ATCGN", "TAGCN")


class GNode:

    """
    Genome Node representing continuous sequence
    """

    def __init__(self, chrom, abs_start, abs_end, is_alt=False, sequence=None, strand=True, complement=None, var_data=None):
        """
        chrom/abs_start/abs_end are the absolute positions of the sequence in the reference
        is_alt is True if the GNode represents variant allele sequence
        If is_alt, the variant sequence must be provided. Otherwise, a biograph.Reference object will be parsed to
                produce sequence.
        strand is True/False for +/-
        complement is a pointer to the Node's opposite strand counterpart
        var_data hold information pertinent to the variant, typically vcf records
        """
        self.chrom = chrom
        self.abs_start = abs_start
        self.abs_end = abs_end
        self.is_alt = is_alt
        if sequence is not None:
            self.sequence = Sequence(sequence) if isinstance(sequence, str) else sequence
        else:
            self.sequence = None
        self.strand = strand
        self.complement = complement
        self.var_data = var_data

        self.cache = dict()

    def get_sequence(self, reference):
        """
        Returns this Node's sequence,
        The reference must be provided if self.sequence is None
        """
        if self.sequence is None:
            # Correction because biograph.Reference is awful...
            return self.__new_ref_worker(reference)

        if self.strand:
            return str(self.sequence)

        return str(self.sequence.rev_comp())

    def get_size(self):
        """
        Get variant's size relative to reference.
        If var_data is a vcf.VCFRecord, leverage that to compute size.
        Otherwise, we'll use the coordinates and
        Deletions have a negative size and insertions have a positive
        """
        if isinstance(self.var_data, vcf.model._Record):  # pylint: disable=protected-access
            #fail = False
            if self.var_data.is_sv:
                ret = None
                try:
                    if isinstance(self.var_data.INFO["SVLEN"], list):
                        ret = self.var_data.INFO["SVLEN"][0]
                    else:
                        ret = self.var_data.INFO["SVLEN"]
                except KeyError:
                    ret = abs(len(self.var_data.REF) - len(self.var_data.ALT))
                return ret
            if self.var_data.var_subtype == 'del':
                return -abs(self.abs_end - self.abs_start)
            return len(self.sequence)
        span = abs(self.abs_end - self.abs_start)
        sq_len = len(self.sequence)
        if sq_len < span:  # del
            return -span

        # SNP/MNP/INS
        return sq_len

    def __new_ref_worker(self, reference, max_cache=1000):
        """
        Turn the biograph.Reference into a continuous sequence
        __ref_worker needs to be refactored into here.
        """
        seqid = "{0}|{1}|{2}".format(self.chrom, self.abs_start, self.abs_end)
        if seqid in self.cache:
            if self.strand:
                return self.cache[seqid]

            return self.cache[seqid].translate(STRREVCMP)[::-1]

        ranges = reference.find_ranges(self.chrom, self.abs_start, self.abs_end)
        ret = ""

        if len(ranges) == 1:
            ret = self.__ref_worker(reference)
        else:
            if ranges[0].start > self.abs_start:
                ret = "N" * (ranges[0].start - self.abs_start)
            i = 0
            while i < len(ranges):
                if ranges[i].end < self.abs_end:
                    ret += str(ranges[i].sequence)
                    # hack - I really need to put some thought into this
                    # Also, I should redo the HG001 all_missing analysis to figure out
                    # what weird call made this happen.
                    if i == len(ranges) - 1:
                        ret += ("N" * (self.abs_end - ranges[i].end))
                    else:
                        ret += ("N" * (ranges[i + 1].start - ranges[i].end))
                else:
                    ret += str(ranges[i].sequence)
                i += 1

        self.cache[seqid] = ret
        # Don't let the cache get too big
        if len(self.cache) > max_cache:
            del self.cache[random.choice(self.cache.keys())]

        if self.strand:
            return ret
        return ret.translate(STRREVCMP)[::-1]

    def __ref_worker(self, reference):
        """
        Simplest case of parsing a single ReferenceRange to create the sequence needed
        """
        ret_seq = reference.make_range(self.chrom, self.abs_start, self.abs_end, use_exact_loci=False)
        if ret_seq.start == self.abs_start and ret_seq.end == self.abs_end:
            ret = str(ret_seq.sequence)
        else:
            # Special handling to make the sequences over Ns
            up_seq = ""
            if ret_seq.start > self.abs_start:
                up_seq = "N" * (ret_seq.start - self.abs_start)
            seq = ret_seq.sequence
            dn_seq = ""
            # doing something unoptimal here
            if ret_seq.start + len(ret_seq.sequence) < self.abs_end:
                dn_seq = "N" * (self.abs_end - (ret_seq.start + len(ret_seq.sequence)))
            ret = up_seq + str(seq) + dn_seq
        return ret

    # Rich comparison needed for bisection
    def __lt__(self, other):
        if self.chrom == other.chrom:
            if self.abs_start == other.abs_start:
                return self.abs_end < other.abs_start
            return self.abs_start < other.abs_start
        return self.chrom < other.chrom

    def __gt__(self, other):
        if self.chrom == other.chrom:
            if self.abs_start == other.abs_start:
                return self.abs_end > other.abs_end

            return self.abs_start > other.abs_start
        return self.chrom > other.chrom

    def __eq__(self, other):
        return self.chrom == other.chrom \
            and self.abs_start == other.abs_start \
            and self.abs_end == other.abs_end

    def __repr__(self):
        return "GNode<%s[%d-%d)%s%s>" % (self.chrom, self.abs_start, self.abs_end,
                                         "alt" if self.is_alt else "ref", "+" if self.strand else "-")


class GenomeGraph:

    """
    Represents a sequence as a graph that can be edited by adding edges for genome variation over alleles

    We begin by creating forward and reverse sequence nodes into the graph for every reference chromosome.
    As variants are added into the graph, nodes are split at the boundaries and edges join the reference allele and new
    variant sequence's node.

    The reference nodes are continually tracked and indexed to allow for efficient fetching of regions of the genome.

    Helper methods exist to allow parsing directly from vcf entries and to query the graph structure to produce
    alternate and reference allele paths through the genome.
    """

    def __init__(self, reference):
        """
        Construct a GenomeGraph - initialize with no variants
        """
        self.graph = nx.DiGraph()
        self.index = []
        self.reference = reference
        self.__load_reference()

    def __load_reference(self):
        """
        Build the reference for all the scaffolds
        """
        scafs = self.reference.scaffold_lens
        for key in scafs:
            self.__add_ref(key, 0, scafs[key])

    def __add_ref(self, chrom, abs_start, abs_end):
        """
        Adds a new section of reference into the graph
        """
        new_node = GNode(chrom, abs_start, abs_end)
        self.graph.add_node(new_node)
        self.index_add(new_node)

        new_node_rev = GNode(chrom, abs_start, abs_end, strand=False, complement=new_node)
        new_node.complement = new_node_rev

        self.graph.add_node(new_node_rev)

    def get_node_idx(self, s_node):
        """
        Return the position in the index of a reference node containing chrom pos
        """
        pos = bisect.bisect_left(self.index, s_node)
        # Correction for between chromosomes
        if 0 < pos < len(self.index) and s_node.chrom != self.index[pos].chrom:
            return pos - 1
        return pos

    def get_node(self, chrom, s_pos):
        """
        Get the node containing this reference position
        """
        s_node = GNode(chrom, s_pos, s_pos)
        pos = self.get_node_idx(s_node)
        if pos == 0:
            return self.index[0]
        if pos >= len(self.index):
            return self.index[-1]

        # Half open thus <= <... why did I do this? because it's bisect left...?
        if self.index[pos - 1].chrom == chrom and self.index[pos - 1].abs_start <= s_pos < self.index[pos - 1].abs_end:
            return self.index[pos - 1]

        return self.index[pos]

    def get_chr_end_node(self, chrom):
        """
        Return the end node of a chromosome
        """
        chr_len = self.reference.scaffold_lens[chrom]
        return self.get_node(chrom, chr_len)

    def index_remove(self, node):
        """
        Remove a node from the index
        returns the position at which this node was found
        """
        pos = bisect.bisect_left(self.index, node)
        if self.index[pos] == node:
            del self.index[pos]
        return pos

    def index_add(self, new_node):
        """
        Insert a node into the index
        """
        bisect.insort_left(self.index, new_node)

    def split_ref_pos(self, chrom, pos):
        """
        At the given location, cut the reference sequence into separate nodes and add an
        edge if it hasn't been split already.
        This is in zero based coordinates, so you're looking between bases when specifying pos
        e.g.
            ATC
        you're adding a variant T>G
            split_ref_pos('chr', 1)
        builds
            A-T-C
        So you can then build an edge between node[0], node[2]

        returns the edge of interest
        """
        logging.debug("Splitting @ %s %d", chrom, pos)
        ref_node = self.get_node(chrom, pos)
        if ref_node.abs_start == pos:
            # edge already exists
            if self.index[0] == ref_node:
                return None, ref_node
            up_node = self.index[self.get_node_idx(ref_node) - 1]
            return up_node, ref_node

        up_node_start = ref_node.abs_start
        up_node_end = pos
        up_node = GNode(chrom, up_node_start, up_node_end)
        up_node_rev = GNode(chrom, up_node_start, up_node_end, strand=False, complement=up_node)
        up_node.complement = up_node_rev

        dn_node_start = pos
        dn_node_end = ref_node.abs_end
        dn_node = GNode(chrom, dn_node_start, dn_node_end)
        dn_node_rev = GNode(chrom, dn_node_start, dn_node_end, strand=False, complement=dn_node)
        dn_node.complement = dn_node_rev

        # connect the cut
        self.graph.add_edge(up_node, dn_node, is_alt=0)
        self.graph.add_edge(dn_node_rev, up_node_rev, is_alt=0)

        # inherit the edges
        for i in self.graph.predecessors(ref_node):
            dat = self.graph.get_edge_data(i, ref_node)
            self.graph.add_edge(i, up_node, **dat)
            self.graph.add_edge(up_node_rev, i.complement, **dat)

        for i in self.graph.successors(ref_node):
            dat = self.graph.get_edge_data(ref_node, i)
            self.graph.add_edge(dn_node, i, **dat)
            self.graph.add_edge(i.complement, dn_node_rev, **dat)

        # trash and index
        self.graph.remove_node(ref_node)
        self.graph.remove_node(ref_node.complement)
        old_pos = self.index_remove(ref_node)
        self.index.insert(old_pos, dn_node)
        self.index.insert(old_pos, up_node)

        return up_node, dn_node

    def add_var(self, chrom, start, end, alt_seq, var_data=None):
        r"""
        Incorporates a variant by splitting a continuous region into three parts.
        upstream -> ref -> dnstream
        ref_sequence[chrom][:start] -> ref_sequence[chrom][start:end] -> ref_sequence[chrom][end:]

        Then a new node with alt_seq is added having edges between upstream and dnstream

        var_data is information that will be kept with the created node

        Example: After creating a new GenomeGraph genG and adding a reference
            >>> genG = GenomeGraph(biograph.Reference("ref_folder"))

        Where ref is a biograph.Reference (or GenomeGraph.FakeReference)
        Were we to add one  the following vars, it'd create these graphs
        where -/\ are edges, continuous nucs are nodes

        #snp
        >>> genG.add_var('chr0', 5, 6, 'T')
        ATCAA-G-CACTA
             \T/

        #del
        >>> genG.add_var('chr0', 5, 6, "")
        ATCAA---CACTA
             \G/

        #ins
        >>> genG.add_var('chr0', 5, 5, "CAT")
        ATCAAG---CACTA
             \CAT/
        """
        # TODO add TLOC support
        if start == end:  # up to one cut
            up_node, dn_node = self.split_ref_pos(chrom, start)
        else:  # up to two cuts
            up_node_s, dn_node_s = self.split_ref_pos(chrom, start)
            logging.debug(dn_node_s)
            # tiny bit wasteful when [start-end) contained in same node
            dn_node_e, dn_node_e = self.split_ref_pos(chrom, end)
            up_node = up_node_s
            dn_node = dn_node_e
            logging.debug("placing variant between %s %s ****", up_node, dn_node)

        re_node_start = start
        re_node_end = end
        re_node = GNode(chrom, re_node_start, re_node_end, is_alt=True, sequence=alt_seq, var_data=var_data)
        re_node_rev = GNode(chrom, re_node_start, re_node_end, is_alt=True, sequence=alt_seq,
                            strand=False, complement=re_node, var_data=var_data)
        re_node.complement = re_node_rev

        # tie it in
        self.graph.add_edge(up_node, re_node, is_alt=1)
        self.graph.add_edge(re_node, dn_node, is_alt=1)

        self.graph.add_edge(dn_node.complement, re_node.complement, is_alt=1)
        self.graph.add_edge(re_node.complement, up_node.complement, is_alt=1)

        # could return the new edges so that the calling method
        # e.g. (add_vcf) could data instead of passing it in?

    def add_vcf(self, vcf_record):
        """
        Parses vcf record and call add_var as needed
        """
        for allele in vcf_record.alleles[1:]:
            if allele.type == "BND":
                self.add_bnd(vcf_record, allele)
                continue
            #WARNING - this is incorrect on vt normalized variants
            if allele.type == "DEL":#FOR <DEL> alts
                if "END" not in vcf_record.INFO:
                    raise KeyError("END not found in vcf entry's INFO")
                end = vcf_record.INFO["END"]
                alt_seq = vcf_record.REF
            else: # this is broken. <DUP> <INV> <INS>
                end = vcf_record.end
                alt_seq = allele.sequence
            self.add_var(vcf_record.CHROM, vcf_record.start, end,
                         alt_seq, var_data=vcf_record)

    def add_bnd(self, vcf_record, allele):
        """
        Parsing for BND events

        WARNING - Intra-chromosomal translocations are not supported and will raise RuntimeErrors
        WARNING - we currently assume all BNDs are balanced. UnBalanced TLOCs (like what you'd find in some cancers),
                are not tested
        """
        start = vcf_record.start
        alt_chr = allele.chr
        alt_pos = allele.pos
        alt_seq = str(allele.connectingSequence)

        if vcf_record.CHROM != alt_chr:
            # Just do all the work here
            # same as what we do for the weird BNDS below
            raise RuntimeError("Alleles with TLOCs are not yet supported (vcf: %s, allele: %s)" %
                               (str(vcf_record), str(allele)))

        # allele orientation remoteOrientation
        # A[5:11[ False True
        # C]4:19] False False
        # ]5:11]C True False
        # [2:86[C True True

        # orientation:
        # The orientation of breakend. If the sequence 3' of the breakend is connected, True, else if the sequence 5' of
        # the breakend is connected, False.

        # remote orientation: (essentially mate information)
        # The orientation of breakend's mate. If the sequence 3' of the breakend's mate is connected, True, else if the
        # sequence 5' of the breakend's mate is connected, False.

        # s   t[p[ piece extending to the right of p is joined after t
        if not allele.orientation and allele.remoteOrientation:
            # no strand switch

            # sorting the positions is a cheat to reuse GenomeGraph.add_var
            # I'm only 90% confident that this won't cause a problem
            # 99% confidence would be a lot of testing
            # 100% confidence would be re-implementing add_var type functionality and making this piece of code do the
            # graph cuts and rejoins based on the variant's definition
            if alt_pos < start:  # tloc... wouldn't matter? ew.
                start, alt_pos = alt_pos, start
            self.add_var(alt_chr, start, alt_pos, alt_seq, var_data=vcf_record)

        # s   t]p] reverse comp piece extending left of p is joined after t
        elif not allele.orientation and not allele.remoteOrientation:
            # cut at pos and alt pos
            # join the up_node from from pos cut with the up_node.complement of alt pos cut
            # up_node of pos_cut joins rev_comp_up_node of alt_cut (alt pos)
            # correction for the alt's position becus ainfiana
            # I hate this
            start += 1
            alt_pos -= 1
            alt_seq = alt_seq[1:]
            pos_cut_up_node, pos_cut_dn_node = self.split_ref_pos(vcf_record.CHROM, start)
            alt_cut_up_node, alt_cut_dn_node = self.split_ref_pos(alt_chr, alt_pos)

            # don't forget the alt sequence that's joining the thing
            re_node_start = start
            re_node_end = alt_pos
            re_node = GNode(alt_chr, re_node_start, re_node_end, is_alt=True, sequence=alt_seq, var_data=vcf_record)
            re_node_rev = GNode(alt_chr, re_node_start, re_node_end, is_alt=True, sequence=alt_seq,
                                strand=False, complement=re_node, var_data=vcf_record)
            re_node.complement = re_node_rev

            self.graph.add_edge(pos_cut_up_node, re_node, data=vcf_record, is_alt=1)
            self.graph.add_edge(re_node, alt_cut_up_node.complement, data=vcf_record, is_alt=1)
            # something...
            # For the complement strand
            self.graph.add_edge(alt_cut_up_node, re_node.complement, data=vcf_record, is_alt=1)
            self.graph.add_edge(re_node.complement, pos_cut_up_node.complement, data=vcf_record, is_alt=1)

        # s   ]p]t piece extending to the left of p is joined before t
        elif allele.orientation and not allele.remoteOrientation:
            # no strand switch
            if alt_pos > start:  # do I care about positioning?
                start, alt_pos = alt_pos, start
            self.add_var(alt_chr, start, alt_pos, alt_seq, var_data=vcf_record)

        # s   [p[t reverse comp piece extending right of p is joined before t
        elif allele.orientation and allele.remoteOrientation:
            # alt_cut_dn_node_complement joins pos_cut_dn_node
            # alt_pos -= 1 #??
            # start -= 1
            alt_seq = alt_seq[1:]
            pos_cut_up_node, pos_cut_dn_node = self.split_ref_pos(vcf_record.CHROM, start)
            alt_cut_up_node, alt_cut_dn_node = self.split_ref_pos(alt_chr, alt_pos)

            # don't forget the alt sequence that's joining the thing
            re_node_start = start
            re_node_end = alt_pos
            re_node = GNode(alt_chr, re_node_start, re_node_end, is_alt=True, sequence=alt_seq, var_data=vcf_record)
            re_node_rev = GNode(alt_chr, re_node_start, re_node_end, is_alt=True, sequence=alt_seq, strand=False,
                                complement=re_node, var_data=vcf_record)
            re_node.complement = re_node_rev

            self.graph.add_edge(alt_cut_dn_node.complement, re_node, data=vcf_record, is_alt=1)
            self.graph.add_edge(re_node, pos_cut_dn_node, data=vcf_record, is_alt=1)

            # For the complement
            self.graph.add_edge(pos_cut_dn_node.complement, re_node.complement, data=vcf_record, is_alt=1)
            self.graph.add_edge(re_node.complement, alt_cut_dn_node, data=vcf_record, is_alt=1)

    def all_paths(self, source, sink, max_paths=None):
        """
        Gives all simple paths from source to sink up to max_paths
        Paths are weighted by is_alt, so the first thing returned has a weight of 0 i.e. the reference path
        """
        return islice(nx.all_simple_paths(self.graph, source, sink), max_paths)

    def ref_path(self, source, sink):
        """
        Gives the path of only the reference
        """
        return islice(self.all_paths(source, sink), 0, 1)

    def alt_paths(self, source, sink, max_paths=None):
        """
        Gives all the paths with alternative allele GNodes
        """
        if max_paths is not None:
            max_paths += 1
        return islice(self.all_paths(source, sink, max_paths), 1, None)

    def get_path_seq(self, path):
        """
        Translates a path of nodes into the sequence it represents
        """
        ret = BytesIO()
        for node in path:
            ret.write(node.get_sequence(self.reference))
        ret.seek(0)
        return ret.read()

    def alt_node_iter(self, chrom=None, start=None, end=None):
        """
        Returns an iterator over all alternate nodes in a graph
        if chrom:start-end is provided, restrict to alt_nodes contained within that region
            (All three must be provided or else the iterator is over all alt nodes)
        Both forward and reverse nodes are returned, so filter based on GNode.strand if you want to keep one strand
        """
        if chrom is not None and start is not None and end is not None:
            src = self.get_node(chrom, start)
            snk = self.get_node(chrom, end)
            for path in self.ref_path(src, snk):
                for node in path[:-1]:  # don't want to go past the end
                    for edge in self.graph.successors(node):
                        if edge.is_alt:
                            yield edge
            # need to traverse the complement to return alt_nodes that end in range
            for path in self.ref_path(snk.complement, src.complement):
                for node in path[:-1]:
                    for edge in self.graph.successors(node):
                        if edge.is_alt:
                            yield edge
        else:
            for i in self.graph.nodes():
                if i.is_alt:
                    yield i

    def alt_edge_iter(self, nodes=None):
        """
        For every node with alternate edges, yield a list of all edges
        This is essentially a collapse of all edges from a ref-nodes with alternate paths

        nodes (iterable container, optional (default= all nodes in self.index)) A container of nodes to iterate

        returns a generator where each yield is as list of tuples of edge info [(src, snk, data on edge), ...]
        """
        if nodes is None:
            nodes = self.index

        for node in nodes:
            has_alt = False
            ret = []
            for snk in self.graph.successors(node):
                data = self.graph.get_edge_data(node, snk)
                if data["is_alt"]:
                    has_alt = True
                ret.append((node, snk, data))
            if has_alt:
                yield ret

    def get_allele_seqs(self, chrom, start, end):
        """
        Generates allele sequences for all paths between reference start and end
        Returns a dictionary with two dictionaries of a set of sequences..
        So for a single variant 'v' in a region, you'd see:
        ret["refs"]['v'] = ['ref_sequence']
        ret["alts"]['v'] = ['alt_sequence']

        Note - this generates all possible paths through the reference. So if you have variants SNPs A and B within
        the window, you'll generate 4 sequences with ref-only, only A, only B, A and B.
        Therefore, the only A sequence will be added to ref["refs"]['B'] and vice-versa for only B
        """
        src = self.get_node(chrom, start)
        snk = self.get_node(chrom, end)
        ret = {"refs": defaultdict(set), "alts": defaultdict(set)}
        # Get the set of all variants in the paths
        all_vars = set()
        for path in self.all_paths(src, snk):
            for i in range(len(path) - 1):
                edge_data = self.graph.get_edge_data(path[i], path[i + 1])
                if edge_data["is_alt"]:
                    all_vars.add(path[i] if path[i].is_alt else path[i + 1])

        # Add sequences for each path to the variant
        # and as ref for every variant not in the path.. still trouble?
        for path in self.all_paths(src, snk):
            cur_vars = set()
            has_ref = False
            for i in range(len(path) - 1):
                edge_data = self.graph.get_edge_data(path[i], path[i + 1])
                if edge_data["is_alt"]:
                    cur_vars.add(path[i] if path[i].is_alt else path[i + 1])
                else:  # non-variant edge in this path
                    has_ref = True
            # English1 I could just fake the src's abs_start and the snk's abs_end?
            seq = self.get_path_seq(path)[start - src.abs_start: end - snk.abs_end]
            # add to ref for variants not in this path (we're supporting its reference allele)
            if has_ref:
                for i in all_vars - cur_vars:
                    ret["refs"][i].add(seq)
            for i in cur_vars:
                ret["alts"][i].add(seq)

        return ret

    def get_node_kmer(self, node, kmer_size=50):
        """
        Given a node - beginning at its 5' end, traverse its paths until sequence
        of kmer_size is produced.
        Returns a list of strings with each path's kmer is produced
        if len(node.sequence) >= kmer_size:
            return node.get_sequence(self.reference)[:kmer_size]
        """
        new_seq = node.get_sequence(self.reference)
        if kmer_size - len(new_seq) <= 0:
            return [new_seq[:kmer_size]]
        ret = []
        for key in self.graph.successors(node):
            for extends in self.get_node_kmer(key, kmer_size - len(new_seq)):
                ret.append(new_seq + extends)
        return ret

    def get_edge_kmer(self, edge, kmer_size=50):
        """
        Given an edge in self.graph, generate a kmer centered around it.
        if kmer_size is an odd number, the returned kmer will be 1bp shorter than kmer_size.
        Note: kmer is always reported in forward orientation
        e.g. graph:
        ATCG-CATG

        >>> edge = (ATCG, CATG)
        >>> get_edge_kmer(edge, 4)
        CGCA
        """
        src, snk = edge
        # If you want the end of a node you get_node_kmer complement and move up-wards
        src_edges = self.get_node_kmer(src.complement, int(kmer_size / 2))
        snk_edges = self.get_node_kmer(snk, int(kmer_size / 2))
        ret = []
        for i in src_edges:
            if i.count("N"):  # Can't use biograph.Sequence
                iadd = i.translate(STRREVCMP)[::-1]
            else:
                iadd = Sequence(i).rev_comp()
            for j in snk_edges:
                ret.append(iadd + j)
        return ret
