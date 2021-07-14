"""
BioGraph Precision Compare

Builds kmer-queries around entries in a vcf and adds supporting evidence read counts per-sample.
Also adds INFO annotation of mendelian consistency
In the future could adds phasing information
    (plus feature would use mates)
    (though that may become another tool all on it's own)

"""
# pylint: disable=invalid-name

from __future__ import print_function

import sys
import json
import argparse
import progressbar

import vcf
import biograph
import biograph.tools.log as log
from biograph.utils import load_regions, Pedigree, _PedSample
from biograph.internal import GenomeGraph

class NodeKmerQuery:

    """
    Object helper for performing a KmerQuery over a GenomeGraph
    This was extracted from PrecisionCompare to make it easier to use in other scripts
    """

    def __init__(self, my_bg, genome):
        self.my_bg = my_bg
        self.genome = genome

    def make_kmers(self, alt_node, kmersize=50):
        r"""
        Given a GNode for an alternate allele,
        Make and return a tuple with the sets of ref and alt kmers spanning the
        upstream and downstream edges of the node and into the variant

        For example, given the insertion of AAT, creating all 4mers would return:

        ATCAAG---CACTA
             \AAT/
        ref: AGCA
        alt: AGAA, ATCA
        """
        def get_kmer_set(edges):
            """
            Given a list of edges, return a set of kmers
            """
            ret = set()
            for edge in edges:
                ret.update(self.genome.get_edge_kmer(edge, kmersize))
            return ret

        # direct alt node - queries variant 3'->5' ref
        dan = list(self.genome.graph.edges(alt_node))
        alt_ret = get_kmer_set(dan)

        # complement alt node - queries ref 3'->5' variant
        can = list(self.genome.graph.edges(alt_node.complement))
        alt_ret.update(get_kmer_set(can))

        # downstream ref node - queries ref under variant 3' -> 5'ref
        drn = [x for x in self.genome.graph.edges(dan[0][1].complement) if not (x[0].is_alt or x[1].is_alt)]
        ref_ret = get_kmer_set(drn)

        # upstream ref node - queries ref 3'->5' ref under variant
        crn = [x for x in self.genome.graph.edges(can[0][1].complement) if not (x[0].is_alt or x[1].is_alt)]
        ref_ret.update(get_kmer_set(crn))
        return ref_ret, alt_ret

    def read_search(self, alt_node, kmersize=50, editdistance=0):
        """
        Given an alternate GNode, return a tuple with sets of ReadmapEntry objects for reads containing the ref, alt
        kmers
        """
        ref_kmers, alt_kmers = self.make_kmers(alt_node, kmersize)
        ref_reads = set()
        for seq in ref_kmers:
            ref_reads.update(self.my_bg.read_search(seq, max_mismatch=editdistance))
        alt_reads = set()
        for seq in alt_kmers:
            alt_reads.update(self.my_bg.read_search(seq, max_mismatch=editdistance))

        return ref_reads, alt_reads

    def read_counts(self, alt_node, kmersize=50, editdistance=0):
        """
        Calls NodeKmerQuery read search and simply returns the Ref and Alt counts
        """
        ref, alt = self.read_search(alt_node, kmersize, editdistance)
        return len(ref), len(alt)


class PrecisionCompare:

    """
    BioGraph Precision Compare tool
    """

    def __init__(self, my_bg, my_vcf, my_ref, my_ped, output):
        self.my_bg = my_bg
        self.my_vcf = my_vcf
        self.my_ref = my_ref
        self.my_ped = my_ped
        self.output = output
        self.genome = GenomeGraph(my_ref)
        self.kmer_query = NodeKmerQuery(self.my_bg, self.genome)

        self.variant_count = 0
        self.pgbar = None

    def load_vcf(self, regions=None):
        """
        Load in all the vcf entries within regions
        Regions are pybedtools Intervals
        """
        if regions is None:
            for entry in self.my_vcf:
                try:
                    self.genome.add_vcf(entry)
                    self.variant_count += 1
                except RuntimeError:
                    log.warning("TLOC variant %s excluded", str(entry))
            log.info("Loaded %d vcf entries", self.variant_count)
        else:
            for region in regions:
                for entry in self.my_vcf.fetch(region.chrom, region.start, region.end):
                    # Fetch gives OVERLAPPING variants, so redundant SVs get put in..
                    # This is one solution, but what if the end is in there..
                    # Should consider just stop using HQRegions... or like pre-filter stuff.
                    # But then the intersect tools do the overlap, too...
                    if not region.start <= entry.start < region.end:
                        continue

                    try:
                        self.genome.add_vcf(entry)
                        self.variant_count += 1
                    except RuntimeError:
                        log.warning("TLOC variant %s excluded", str(entry))
            log.info("Loaded %d vcf entries", self.variant_count)


    def annotate(self, kmersize=50, min_size=0, editdistance=0, ref_k_matches=False, kmer_seq_out=False):
        """
        Annotate every entry in the vcf
        kmersize - size of the kmer query to create
        editdistance (not implemented) - maximum editdistance between kmer queries and reads
        alt_refmatch - maximum number of reference matches the alt kmers can have
        max_refmatch - maximum number of reference matches the ref kmers can have
        min_size (not implemented) - minimum size of variants to annotate (min_size=2 excludes SNPs)
        """

        # TODO: Doing this per sample in the pedigree (readmap switching and all that)
        # TODO: when you have neighboring variants, if one is real and the other is fake the
        #      real count needs to be added to PR
        # TODO: We can add phasing information
        # TODO: Get edit_distance kmer query working and replace my_bg.find (masking is important)

        processed = 0
        for alt_node in self.genome.alt_node_iter():
            if not alt_node.strand:
                continue

            if abs(alt_node.get_size()) < min_size:
                # self.output.write_record(vcf_entry) Might should be doing this
                continue

            vcf_entry = alt_node.var_data
            self.genotype_format_prep(vcf_entry)
            gt_lookup = {}
            n_call_data = None
            # I'm going to need this thingy
            for v_sample in vcf_entry.samples:
                if n_call_data is None:
                    try:
                        # If a site has two alleles, they make two nodes
                        # we're writing the nodes redundantly
                        n_call_data = vcf.model.make_calldata_tuple(
                            v_sample.data._fields + ("PR", "PA", "PGT", "PGQ"))
                    except ValueError:
                        continue
                if v_sample.sample not in self.my_ped:
                    v_sample.data += ('.', '.', '.', '.')
                    v_sample.data = n_call_data(*v_sample.data)
                    continue

                n_data = self.genotype_anno(v_sample, alt_node, kmersize, editdistance, n_call_data)
                #n_data = None;

                if n_data is not None:
                    gt_lookup[v_sample.sample] = n_data

            if ref_k_matches or kmer_seq_out:
                self.kmer_sequence_anno(vcf_entry, alt_node, kmersize, ref_k_matches, kmer_seq_out)
            # de-novo annotation Warning! this doesn't actually work. Going to need to post-process it
            # self.inheritance_anno(vcf_entry, gt_lookup)

            self.output.write_record(vcf_entry)
            processed += 1
            if not processed % 10:
                self.pgbar.update(min(processed, self.variant_count))

    def genotype_format_prep(self, vcf_entry):  # pylint: disable=no-self-use
        """
        Adds the relevant FORMAT information to a sample
        """
        # Need to get the fields correct
        # there's something weird happening where this is called twice somehow sometimes or something
        if vcf_entry.FORMAT is None:
            vcf_entry.FORMAT = "PR"
        fmt_data = vcf_entry.FORMAT.split(':')
        if "PR" not in fmt_data:
            vcf_entry.add_format("PR")
        if "PA" not in fmt_data:
            vcf_entry.add_format("PA")
        if "PGT" not in fmt_data:
            vcf_entry.add_format("PGT")
        if "PGQ" not in fmt_data:
            vcf_entry.add_format("PGQ")

    def genotype_anno(self, v_sample, alt_node, kmersize, editdistance, n_call_data):
        """
        Annotates a sample in a vcf entry with
        """
        self.my_bg.set_readmap(v_sample.sample)
        # TODO: handling compound het...
        # pylint: disable=broad-except
        try:
            ref_count, alt_count = self.kmer_query.read_counts(alt_node, kmersize, editdistance)
        except Exception:
            return None

        gt, gq = biograph.genotyper(ref_count + alt_count, alt_count)
        n_data = (str(ref_count), str(alt_count), gt, '%.2f' % gq)
        v_sample.data += n_data
        v_sample.data = n_call_data(*v_sample.data)
        return n_data

    def kmer_sequence_anno(self, vcf_entry, alt_node, kmersize, ref_k_matches, kmer_seq_out):
        """
        In-place edits the vcf_entry by adding kmer sequences and reference-matching counts of kmers
        """
        try:
            ref_seqs, alt_seqs = self.kmer_query.make_kmers(alt_node, kmersize)
        except RuntimeError as e:
            log.debug('could not add anno on %s', str(vcf_entry))
            log.debug(e)
            return
        if ref_k_matches:
            if "PRM" not in vcf_entry.INFO:
                vcf_entry.add_info("PRM", [])
            if "PAM" not in vcf_entry.INFO:
                vcf_entry.add_info("PAM", [])
            ref_cnt = 0
            alt_cnt = 0
            for i in ref_seqs:
                i = biograph.Sequence(i) if isinstance(i, biograph.Sequence) else i
                ref_cnt = self.my_ref.find(i).matches + self.my_ref.find(i.rev_comp()).matches
                vcf_entry.INFO['PRM'].append(ref_cnt)
            for i in alt_seqs:
                i = biograph.Sequence(i) if isinstance(i, biograph.Sequence) else i
                alt_cnt = self.my_ref.find(i).matches + self.my_ref.find(i.rev_comp()).matches
                vcf_entry.INFO['PAM'].append(alt_cnt)
        if kmer_seq_out:
            if "PRefK" not in vcf_entry.INFO:
                vcf_entry.add_info("PRefK", [])
            if "PAltK" not in vcf_entry.INFO:
                vcf_entry.add_info("PAltK", [])
            for i in ref_seqs:
                vcf_entry.INFO["PRefK"].append(i)
            for i in alt_seqs:
                vcf_entry.INFO["PAltK"].append(i)

    def inheritance_anno(self, vcf_entry, gt_lookup):
        """
        In-place edits the vcf_entry by adding the Inheritance annotation
        """
        for sample in self.my_ped:
            sample = self.my_ped[sample]
            # Has parents, so we can potentially work it in
            if sample.mother and sample.father:
                if sample.ind_id not in gt_lookup:  # Wasn't genotyped pass?
                    vcf_entry.add_info("Inheritance", "NeedsSquared")
                elif gt_lookup[sample.ind_id][2] == '0/0':
                    if gt_lookup[sample.mother.ind_id][0] == '1/1' and gt_lookup[sample.father.ind_id][2] == '1/1':
                        vcf_entry.add_info("Inheritance", "Inconsistent")
                elif gt_lookup[sample.ind_id][2] == '0/1':
                    if gt_lookup[sample.mother.ind_id][2] == '1/1' and gt_lookup[sample.father.ind_id][2] == '1/1':
                        vcf_entry.add_info("Inheritance", "Inconsistent")
                    elif gt_lookup[sample.mother.ind_id][2].endswith('/1') or gt_lookup[sample.father.ind_id][2].endswith('/1'):
                        pass
                    else:
                        vcf_entry.add_info("Inheritance", "DeNovoCandidate")
                elif gt_lookup[sample.ind_id][2] == '1/1':
                    if gt_lookup[sample.mother.ind_id][2] == '0/0' and gt_lookup[sample.father.ind_id][2] == '0/0':
                        vcf_entry.add_info("Inheritance", "DeNovoCandidate")
                # Need . alleles
                #
                # vcf_entry.add_info("Inheritance", "CoverageIssue")

    # def __hide_map_filter(self, ctx, threshold, alt_ctx, ref_ctx):
        # Forget the mapping stuff, for now
        # Will count identical reads as well as sub-reads
        # We could do a de-duping here if we wanted
        # Wait... the matches filtering.... also
        # def map_filter(ctx, threshold):
            # need to do stranded-ness...
            # return not (self.my_ref.find(ctx.sequence).matches > threshold or
                        # self.my_ref.find(ctx.sequence.rev_comp()).matches > threshold)
        # and should we be doing an exclusion? like, if there's the same entry in alt and ref?
        # should make this optional - like if it's off (-1) get out of there
        # alt_depth = sum([len(self.my_bg.entry_to_reads(x))
                         # for x in alt_ctx - ref_ctx
                         # if map_filter(x, alt_refmatch)])
        # ref_depth = sum([len(self.my_bg.entry_to_reads(x))
                         # for x in ref_ctx - alt_ctx
                         # if map_filter(x, max_refmatch)])

        # return ref_depth, alt_depth


def parse_args(args):
    """ Make pretty arguments """
    parser = argparse.ArgumentParser(description='Run kmer analysis with variants from proband')
    # io group
    ioGroup = parser.add_argument_group("I/O Arguments")
    ioGroup.add_argument("-b", "--biograph", metavar="BG", required=True,
                         help="Merged BioGraph file containing individuals")
    ioGroup.add_argument("-v", "--variant_file", metavar="VCF", required=True,
                         help="The Anchored Assembly vcf containing all samples")
    ioGroup.add_argument("-r", "--ref", metavar="REF", required=True,
                         help="Reference genome folder")
    ioGroup.add_argument("-p", "--ped", metavar="PED", required=False,
                         help="Pedigree file with individuals to be processed")
    ioGroup.add_argument("-o", "--output", metavar="OUT", default="output.vcf",
                         help="Output annotated vcf file.")
    # Tweaking params
    twGroup = parser.add_argument_group("Analysis Arguments")
    twGroup.add_argument("-k", "--kmersize", default=50, type=int,
                         help="Number of basepairs used in the kmer query (%(default)s)")
    twGroup.add_argument("-e", "--editdistance", default=0, type=int,
                         help="Maximum mismatches between a read and the kmer used (%(default)s)")
    # twGroup.add_argument("-a", "--alt-refmatch", default=0, type=int,
                         # help="Maximum number of reference matches for alt kmer queries (%(default)s)")
    # twGroup.add_argument("-m", "--max-refmatch", default=1, type=int,
                         # help="Maximum number of reference matches for ref kmer queries (%(default)s)")

    # Filtering params - what to include/exclude from the analysis
    flGroup = parser.add_argument_group("Filtering Arguments")
    flGroup.add_argument("-s", "--size-min", default=1, type=int,
                         help="Minimum span of variants to add support. 2 would exclude SNPs (%(default)s)")
    flGroup.add_argument("-R", "--region_bed", metavar="BED", required=False,
                         help=("Bed file with regions within which to analyze variants. Only events within POS "
                               "in regions with be included"))
    flGroup.add_argument("--region", metavar="REG", required=False, nargs="*",
                         help=("Region over which to analyze variants in form chrom:start-end. Only events "
                               "within POS in regions with be included"))
    flGroup.add_argument("-i", "--individual", metavar="ID", nargs="+", required=False,
                         help="Restrict analysis to individual ids")
    flGroup.add_argument("-I", "--Individual", metavar="ID", nargs="+", required=False,
                         help="Exclude individuals from analysis")
    flGroup.add_argument("-f", "--family", metavar="ID", nargs="+", required=False,
                         help="Restrict analysis to family ids")
    flGroup.add_argument("-F", "--Family", metavar="ID", nargs="+", required=False,
                         help="Exclude family ids from analysis")

    anGroup = parser.add_argument_group("Annotation Arguments")
    anGroup.add_argument("--ref-counts", action="store_true",
                         help="Add the PRM, PAM for kmer's reference matching counts")
    anGroup.add_argument("--kmer-seq", action="store_true",
                         help="Write Kmers used to query reads")
    rnGroup = parser.add_argument_group("Runtime Arguments")
    rnGroup.add_argument("--debug", action="store_true",
                         help="Verbose logging")
    args = parser.parse_args(args)
    log.setup_logging(args.debug)

    log.debug("Params:\n%s", json.dumps(vars(args), indent=4))

    return args


def edit_vcf_header(vcf_file, ref_counts=False, kmer_seq=False):
    """
    In-place addition of the new Formats and Infos into the vcf_file
    """
    # pylint: disable=protected-access
    vcf_file.formats['PR'] = vcf.parser._Format(
        id="PR", num=1, type='Integer', desc="PCmp Reference Allele Supporting Reads")
    vcf_file.formats['PA'] = vcf.parser._Format(
        id="PA", num=1, type='Integer', desc="PCmp Alternate Allele Supporting Reads")
    vcf_file.formats['PGT'] = vcf.parser._Format(id="PGT", num=1, type='String', desc="PCmp Genotype")
    vcf_file.formats['PGQ'] = vcf.parser._Format(id="PGQ", num=1, type='Float', desc="PCmp Genotype Quality")
    vcf_file.infos['Inheritance'] = vcf.parser._Info(
        id="Inheritance", num=1, type='String', desc="Annotation on inheritance pattern observed", source=None,
        version=None)
    if ref_counts:
        vcf_file.infos["PRM"] = vcf.parser._Info(source=None, version=None, num=".",
                                                 id="PRM", type='Integer', desc="PCmp ref-kmer matches in reference")
        vcf_file.infos["PAM"] = vcf.parser._Info(source=None, version=None, num=".",
                                                 id="PAM", type='Integer', desc="PCmp alt-kmer matches in reference")
    if kmer_seq:
        vcf_file.infos["PRefK"] = vcf.parser._Info(source=None, version=None, num=".",
                                                   id="PRefK", type='String', desc="PCmp ref-kmer sequences used in query")
        vcf_file.infos["PAltK"] = vcf.parser._Info(source=None, version=None, num=".",
                                                   id="PAltK", type='String', desc="PCmp alt-kmer sequences used in query")


def kpcmp_run(args):
    """
    Parse args and call PCmp
    """
    args = parse_args(args)
    my_vcf = vcf.Reader(open(args.variant_file))
    my_ref = biograph.Reference(args.ref)
    my_ped = None
    my_bg = biograph.BioGraph(args.biograph, biograph.CacheStrategy.MMAPCACHE)
    if args.ped is not None:
        my_ped = Pedigree(args.ped)
        my_ped.filter(args.family, args.Family, args.individual, args.Individual)
    else:
        if len(my_bg.metadata.samples) != 1 or len(my_vcf.samples) != 1:
            raise RuntimeError("A pedigree file must be provided for BioGraphs or VCFs with more than one sample")
        sample_name = my_vcf.samples[0]
        # This will break if we do advanced stuff on the Pedigree object
        my_ped = {sample_name: _PedSample("0", sample_name, "0", "0", "0", "0")}

    edit_vcf_header(my_vcf, args.ref_counts, args.kmer_seq)
    my_out = vcf.Writer(open(args.output, 'w'), my_vcf)

    pcmp = PrecisionCompare(my_bg, my_vcf, my_ref, my_ped, my_out)
    regions = load_regions(args.region_bed, args.region)

    # Load in the entries
    log.info("Loading VCF Entries")
    pcmp.load_vcf(regions)
    # Output the entries
    log.info("Annotating Entries")

    pcmp.pgbar = progressbar.ProgressBar(redirect_stdout=True, max_value=pcmp.variant_count, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(pcmp.variant_count), '] ',
        progressbar.Bar(), ' (', progressbar.ETA(), ') ',
    ])

    pcmp.annotate(args.kmersize, args.size_min, args.editdistance, args.ref_counts, args.kmer_seq)
    pcmp.pgbar.finish()

    log.info("Finished")


if __name__ == "__main__":
    kpcmp_run(sys.argv[1:])
