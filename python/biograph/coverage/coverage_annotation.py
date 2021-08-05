"""
Calculates read coverage attributes for VCF entries
TODO: Unite with precision_comparison in the same style as GTAnno
"""
import logging

from collections import OrderedDict
from statistics import median
import itertools

import numpy as np

import biograph.variants as bgexvar
from biograph.coverage import make_annotation_base

AlleleCoverageBase = make_annotation_base("AlleleCoverageBase", [
    # These are based on breakpoints
    ("upstream_span", "US", "Length of best-spanning read upstream", 0),
    ("dnstream_span", "DS", "Length of best-spanning read dnstream", 0),
    ("upstream_coverage", "UC", "Number of reads spanning upstream", 0),
    ("dnstream_coverage", "DC", "Number of reads spanning dnstream", 0),
    ("upstream_dir_cov", "UDC", "Upstream direct-strand coverage", 0),
    ("upstream_com_cov", "UCC", "Upstream complement-strand coverage", 0),
    ("dnstream_dir_cov", "DDC", "Dnstream direct-strand coverage", 0),
    ("dnstream_com_cov", "DCC", "Dnstream complement-strand coverage", 0),
    ("upstream_min_overlap", "UMO", "Upstream minimum read-overlap", np.iinfo(np.int16).max),
    ("dnstream_min_overlap", "DMO", "Dnstream minimum read-overlap", np.iinfo(np.int16).max),
    ("upstream_max_overlap", "UXO", "Upstream maximum read-overlap", 0),
    ("dnstream_max_overlap", "DXO", "Dnstream maximum read-overlap", 0),

    # These are over entire assemblies
    ("num_reads", "NR", "Number of reads supporting assembly", 0),
    ("min_overlap", "MO", "Minimum overlap of reads supporting assembly", np.iinfo(np.int16).max),
    ("max_overlap", "XO", "Maximum overlap of reads supporting assembly", 0),
    ("max_coverage", "XC", "Maximum coverage over assembly", 0),
    ("average_coverage", "AC", "Average coverage over assembly", 0),
    ("median_coverage", "EC", "Median coverage over assembly", 0),
    ("minimum_coverage", "MC", "Minimum coverage over assembly", 0),
    ("read_cov_max_paths", "MP", "coverage_max_paths needed to get all the coverage on this assembly", 0),

    # Total from align counts
    ("local_read_lens", "AC_LR", "Align count sum of distinct reads of local read lengths", 0),
    ("local_aligned_bases", "AC_LA", "Align count sum of local read lengths", 0),
    ("tot_aligned_bases", "AC_TA", "Align count sum of operlapping aligned lengths", 0)
])

# Disable some warnings that aren't compatible with make_annotation_base.
# pylint: disable=attribute-defined-outside-init,too-many-instance-attributes,assigning-non-slot

class AlleleCoverage(AlleleCoverageBase):
    """
    Holds all the coverage information for a single allele
    """

    __slots__ = ["asm", "left_offset", "right_offset", "matches_reference", "assembly_id", "read_map"]

    def __init__(self, asm, read_map):
        """
        Hold attibutes
        """
        super().__init__()
        self.asm = asm
        self.left_offset = asm.left_offset
        self.right_offset = asm.right_offset
        self.matches_reference = asm.matches_reference
        self.assembly_id = asm.assembly_id
        self.read_map = read_map

    def parse_reads(self):
        """
        Intake the reads from a particular assembly
        """
        cov = self.asm.read_coverage
        seq_len = len(self.asm.seq)
        # interbase coverage
        interbase_depth = cov.calc_depths()
        interbase_dir = cov.calc_depths(include_rev=False, readmap=self.read_map)
        self.add_from_depths(seq_len, cov, interbase_depth, interbase_dir)

    def add_ref_depths(self, ref_is_upstream, cov, interbase_depth, interbase_dir):
        """
        Save information on reference coverage to all variants that are ref_is_upstream
        """
        self.num_reads += cov.get_tot_read_count()

        # Calculate coverage depths
        self.max_coverage = max(interbase_depth)
        self.average_coverage = round(sum(interbase_depth)/len(interbase_depth))
        self.median_coverage = round(median(interbase_depth))
        self.minimum_coverage = min(interbase_depth)

        if ref_is_upstream:
            self.upstream_coverage = interbase_depth[0]
            self.upstream_dir_cov = interbase_dir[0]
            self.upstream_com_cov = interbase_depth[0] - interbase_dir[0]
        else:
            self.dnstream_coverage = interbase_depth[0]
            self.dnstream_dir_cov = interbase_dir[0]
            self.dnstream_com_cov = interbase_depth[0] - interbase_dir[0]

        all_overlaps = cov.get_overlap_min_max()
        self.min_overlap = all_overlaps[0]
        self.max_overlap = all_overlaps[1]

        if ref_is_upstream:
            upstream_cov = cov.get_reads_spanning_offset(0)
            self.upstream_span = upstream_cov.get_max_flank(0)
            upstream_overlaps = upstream_cov.get_overlap_min_max()
            self.upstream_min_overlap = upstream_overlaps[0]
            self.upstream_max_overlap = upstream_overlaps[1]
        else:
            dnstream_cov = cov.get_reads_spanning_offset(0)
            self.dnstream_span = dnstream_cov.get_max_flank(0)
            dnstream_overlaps = dnstream_cov.get_overlap_min_max()
            self.dnstream_min_overlap = dnstream_overlaps[0]
            self.dnstream_max_overlap = dnstream_overlaps[1]

    def add_var_depths(self, seq_len, cov, interbase_depth, interbase_dir):
        """
        # Calculate everything we can about this variant without
        # knowing anything about assemblies around it.
        """
        self.num_reads += cov.get_tot_read_count()

        # Calculate coverage depths
        self.max_coverage = max(interbase_depth)
        self.average_coverage = round(sum(interbase_depth)/len(interbase_depth))
        self.median_coverage = round(median(interbase_depth))
        self.minimum_coverage = min(interbase_depth)

        self.upstream_coverage = interbase_depth[0]
        self.upstream_dir_cov = interbase_dir[0]
        self.upstream_com_cov = interbase_depth[0] - interbase_dir[0]

        self.dnstream_coverage = interbase_depth[-1]
        self.dnstream_dir_cov = interbase_dir[-1]
        self.dnstream_com_cov = interbase_depth[-1] - interbase_dir[-1]

        all_overlaps = cov.get_overlap_min_max()
        self.min_overlap = all_overlaps[0]
        self.max_overlap = all_overlaps[1]

        upstream_cov = cov.get_reads_spanning_offset(0)
        self.upstream_span = upstream_cov.get_max_flank(0)
        upstream_overlaps = upstream_cov.get_overlap_min_max()
        self.upstream_min_overlap = upstream_overlaps[0]
        self.upstream_max_overlap = upstream_overlaps[1]

        dnstream_cov = cov.get_reads_spanning_offset(seq_len)
        self.dnstream_span = dnstream_cov.get_max_flank(seq_len)
        dnstream_overlaps = dnstream_cov.get_overlap_min_max()
        self.dnstream_min_overlap = dnstream_overlaps[0]
        self.dnstream_max_overlap = dnstream_overlaps[1]

    def join_asm_attributes(self, other):
        """
        Inherit attributes between self and other. Generaly the max of the two
        """
        self.num_reads = max(self.num_reads, other.num_reads)
        self.min_overlap = min(self.min_overlap, other.min_overlap)
        self.max_overlap = max(self.max_overlap, other.max_overlap)
        self.max_coverage = max(self.max_coverage, other.max_coverage)
        self.average_coverage = max(self.average_coverage, other.average_coverage)
        self.median_coverage = max(self.average_coverage, other.median_coverage)
        self.minimum_coverage = max(self.average_coverage, other.median_coverage)

    def left_join_upstream(self, other):
        """
        Inherit the other's upstream stats (and the unique ones)
        """
        self.upstream_span = other.upstream_span
        self.upstream_coverage = other.upstream_coverage
        self.upstream_dir_cov = other.upstream_dir_cov
        self.upstream_com_cov = other.upstream_com_cov
        self.upstream_min_overlap = other.upstream_min_overlap
        self.upstream_max_overlap = other.upstream_max_overlap

    def left_join_dnstream(self, other):
        """
        Inherit the other's dnstream stats
        """
        self.dnstream_span = other.dnstream_span
        self.dnstream_coverage = other.dnstream_coverage
        self.dnstream_dir_cov = other.dnstream_dir_cov
        self.dnstream_com_cov = other.dnstream_com_cov
        self.dnstream_min_overlap = other.dnstream_min_overlap
        self.dnstream_max_overlap = other.dnstream_max_overlap

    def __str__(self):
        return ("asm_id:{0} is_ref:{9} left_pos:{7} right_pos:{8} upstream_span:{1} "
                "dnstream_span:{2} upstream_coverage:{3} dnstream_coverage:{4} "
                "upstream_stranded_coverage:{5} dnstream_stranded_coverage:{6} num_reads:{11} "
                "min_overlap:{12} max_overlap:{13} max_coverage:{14} median_coverage:{15} "
                "average_coverage:{10} upstream_min_overlap:{16} upstream_max_overlap:{17} "
                "dnstream_min_overlap:{18} dnstream_max_overlap:{19}"
                ).format(self.assembly_id, self.upstream_span,
                         self.dnstream_span, self.upstream_coverage,
                         self.dnstream_coverage,
                         str((self.upstream_dir_cov, self.upstream_com_cov)),
                         str((self.dnstream_dir_cov, self.dnstream_com_cov)),
                         self.left_offset, self.right_offset,
                         self.matches_reference, self.average_coverage,
                         self.num_reads, self.min_overlap, self.max_overlap,
                         self.max_coverage, self.median_coverage,
                         self.upstream_min_overlap, self.upstream_max_overlap,
                         self.dnstream_min_overlap, self.dnstream_max_overlap)


class CovAnno:
    """
    Calculate the Coverage Annotation SAMPLE fields, can return a dict
    that will be added to a vcf line

    We can make a new *Anno object as long as it has
    1) def parse(self, rmap, ref_file, chrom, entries):
    2) staticmethod - get_header(): list of strings of vcf header entries
    3) staticmethod - get_format_tags(): list of FORMAT column tags
    4) staticmethod - get_blank_format(): returns the dict of default values generated for f
    """

    def __init__(self, min_insert=200, max_insert=1000, max_reads_per_entry=0, max_coverage_paths=0, phasing=False, max_coverage_alleles=None, filter_dup_align=False, ideal_insert=None, placer_max_ambig=15):
        """
        Coverage properties over alleles
        """
        self.min_dist = min_insert
        self.max_dist = max_insert
        self.max_reads_per_entry = max_reads_per_entry
        self.max_coverage_paths = max_coverage_paths
        self.readmap = None
        self.phasing = phasing
        self.max_coverage_alleles = max_coverage_alleles
        self.filter_dup_align = filter_dup_align
        self.ideal_insert = ideal_insert
        self.placer_max_ambig = placer_max_ambig

    def ref_alt_locus_uniter(self, read_map, asms): # pylint: disable=too-many-statements,no-self-use
        """
        pairs alt assemblies with their correct ref_alt
        yields tuple of (ref AlleleCoverage, alt AlleleCoverage)
        """
        def cleanup(lookup, position):
            """
            clear out the entries in lookup ending before position
            """
            keep_popping = True
            while keep_popping:
                keep_popping = False
                m_pos_l, m_pos_r = next(iter(lookup[0])), next(iter(lookup[1]))
                if m_pos_l and m_pos_l < position:
                    del(lookup[0][m_pos_l])
                    keep_popping = True
                if m_pos_r and m_pos_r < position:
                    del(lookup[1][m_pos_r])
                    keep_popping = True

        def make_ref_cov(ref_lookup, alt):
            """
            Make a corresponding ref cov for this alt
            """
            left_ref = ref_lookup[0][alt.asm.left_offset]
            ref_allele = AlleleCoverage(left_ref.asm, None)  # Make a blank one
            ref_allele.left_offset = left_ref.asm.left_offset
            ref_allele.read_cov_max_paths = left_ref.asm.read_cov_max_paths
            try:
                right_ref = ref_lookup[1][alt.asm.right_offset]
            except KeyError:
                logging.error("Forcing left_ref @ position %d for %s.", alt.asm.right_offset, alt)
                return left_ref

            ref_allele.read_cov_max_paths = max(ref_allele.read_cov_max_paths, right_ref.asm.read_cov_max_paths)

            # point insertions are weird
            if alt.asm.left_offset == alt.asm.right_offset:
                left_ref, right_ref = right_ref, left_ref
                ref_allele.left_offset = right_ref.left_offset
                ref_allele.right_offset = left_ref.right_offset
                # Just take whichever breakpoint has more coverage
                # point information.. should be identical to downstream
                cov_fil = max(left_ref.upstream_coverage, right_ref.dnstream_coverage)
                ref_allele.num_reads = cov_fil
                ref_allele.max_coverage = cov_fil
                ref_allele.median_coverage = cov_fil
                ref_allele.average_coverage = cov_fil
                ref_allele.minimum_coverage = cov_fil
                ref_allele.min_overlap = min(left_ref.upstream_min_overlap, right_ref.dnstream_min_overlap)
                ref_allele.max_overlap = max(left_ref.upstream_max_overlap, right_ref.dnstream_max_overlap)
            else:
                ref_allele.left_offset = left_ref.left_offset
                ref_allele.right_offset = right_ref.right_offset
                # Need this MARKERMAK -- Do I want to get the MAX for left and right... especially if they're different...?
                # Also, if the left/right is the same, do I really need to go through all this work?...
                ref_allele.join_asm_attributes(right_ref)

            ref_allele.left_join_upstream(left_ref)
            ref_allele.left_join_dnstream(right_ref)
            return ref_allele

        alt_lookup = []
        ref_lookup = [OrderedDict(), OrderedDict()]  # I want the left and right up ref numbers
        ref_position = 0  # How far down the tracer we've reached
        cur_ref = []
        for asm in asms:
            alt_cov = AlleleCoverage(asm, read_map)
            if asm.matches_reference:
                ref_lookup[0][asm.left_offset] = alt_cov
                ref_lookup[1][asm.right_offset] = alt_cov
                # Parse alts once we're past their ref pos
                ref_position = asm.left_offset - 1
                # Save the reference assembly to feed to subsequent stages.
                cur_ref.append(alt_cov)
                continue
            alt_lookup.append(alt_cov)
            # Alt is ready?
            clear_clean = None
            while alt_lookup and alt_lookup[0].asm.right_offset < ref_position:
                alt_allele = alt_lookup.pop(0)
                if clear_clean is None:
                    # clean everything before the first parsed position
                    clear_clean = alt_allele.left_offset - 1
                ref_cov = make_ref_cov(ref_lookup, alt_allele)
                while cur_ref and cur_ref[0].asm.left_offset < ref_cov.left_offset:
                    yield cur_ref[0], None
                    cur_ref = cur_ref[1:]
                yield ref_cov, alt_allele

            # Clean up reference lookups that end before ref_position
            if clear_clean is not None:
                cleanup(ref_lookup, clear_clean)

        for alt_allele in alt_lookup:
            ref_cov = make_ref_cov(ref_lookup, alt_allele)
            while cur_ref and cur_ref[0].asm.left_offset < ref_cov.left_offset:
                yield cur_ref[0], None
                cur_ref = cur_ref[1:]
            yield ref_cov, alt_allele
        for ref_allele in cur_ref:
            yield ref_allele, None

    def build_vcf(self, ref, alt):  # pylint: disable=no-self-use
        """
        Collapses the ref and alt annotations, populates the info into
        ref.vcf_entry_info.vcf_entry.new_fmt
        yields a reconstructed bex
        """
        # modify the alt first alt - only one with vcf_entry_info
        # TODO... make the pairs for these that are all Number='R'
        combo = {}
        ref_d = ref.to_samp_dict()
        alt_d = alt.to_samp_dict()

        for key in ref_d:
            if isinstance(alt_d[key], float):
                ref_d[key] = "%.2f" % (ref_d[key])
                alt_d[key] = "%.2f" % (alt_d[key])
            combo[key] = "{0},{1}".format(str(ref_d[key]), str(alt_d[key]))
        alt.asm.vcf_entry_info.new_fmt.update(combo)
        return alt.asm

    def with_depths(self, asm):
        """
        Calculate depths of assemblies
        """
        asm.read_depths = asm.read_coverage.calc_depths()
        asm.fwd_read_depths = asm.read_coverage.calc_depths(include_rev=False, readmap=self.readmap)
        asm.allele_coverage = AlleleCoverage(asm, self.readmap)
        asm.ref_allele_coverage = AlleleCoverage(asm, self.readmap)
        asm.allele_coverage.read_cov_max_paths = asm.read_cov_max_paths

        # Calculate everything we can about this variant without
        # knowing anything about assemblies around it.
        if not asm.matches_reference:
            asm.allele_coverage.add_var_depths(len(asm.seq), asm.read_coverage,
                                               asm.read_depths, asm.fwd_read_depths)
        return asm

    def annotate_from_edge(self, ref_offset, left_asms, inserts, right_asms):
        """annotate_from_edge is called for each junction between assemblies.

        left_asms is a list of all the assemblies upstream of this junction,
        right_asms is a list of all the assemblies downstream of this unction,
        and inserts is a list of all the point inserts at this junction
        """

        # Separate out reference assemblies from variant assemblies.
        left_refs = [asm for asm in left_asms if asm.matches_reference]
        right_refs = [asm for asm in right_asms if asm.matches_reference]
        left_vars = [asm for asm in left_asms if not asm.matches_reference]
        right_vars = [asm for asm in right_asms if not asm.matches_reference]

        left_ref_coverage = None
        # Gather coverage from upstream reference assemblies
        for ref_asm in left_refs:
            if left_ref_coverage:
                logging.error(f"Duplicate reference assemblies encountered at {ref_offset}")
            left_ref_coverage = ref_asm.read_coverage.get_and_adjust_reads_spanning_offset(len(ref_asm.seq))
        # Gather coverage from downstream reference assemblies
        right_ref_coverage = None
        for ref_asm in right_refs:
            if right_ref_coverage:
                logging.error(f"Duplicate reference assemblies encountered at {ref_offset}")
            right_ref_coverage = ref_asm.read_coverage.get_and_adjust_reads_spanning_offset(0)

        if left_ref_coverage is None or right_ref_coverage is None:
            if left_vars or right_vars or inserts:
                logging.error(f"missing reference coverage at {ref_offset}")
            return

        # Only count reference reads that continue both upstream *and*
        # downstream in reference.
        ref_coverage = left_ref_coverage.intersection_with(right_ref_coverage)

        # ref_coverage now contains all reads that cross this
        # reference position and continue into reference in both
        # directions.  These reads are cenetered around "0".

        ref_depths = ref_coverage.calc_depths()
        fwd_ref_depths = ref_coverage.calc_depths(include_rev=False, readmap=self.readmap)
        if len(ref_depths) != 1:
            logging.error("We should only have coverage for the interbase position at the edge of the assemblies")
            return

        # Save information on reference coverage to all variants that are upstream.
        for asm in itertools.chain(left_vars, inserts):
            asm.ref_allele_coverage.add_ref_depths(False, # Ref is not upstream
                                                   ref_coverage, ref_depths, fwd_ref_depths)

        # Save information on reference coverage to all variants that are downstream
        for asm in itertools.chain(inserts, right_vars):
            asm.ref_allele_coverage.add_ref_depths(True, # Ref is to the left of asm
                                                   ref_coverage, ref_depths, fwd_ref_depths)

        asm.ref_allele_coverage.read_cov_max_paths = max((a.read_cov_max_paths for a in itertools.chain(left_refs, right_refs)), default=0)

    @staticmethod
    def extract_cov_pass_ref(entries):
        """
        Yields tuples of reference, alternate (if exists) coverages
        """
        for asm in entries:
            var_cov = asm.allele_coverage
            asm.allele_coverage = None
            ref_cov = asm.ref_allele_coverage
            asm.ref_allele_coverage = None

            if asm.matches_reference:
                yield (asm, None)
            else:
                yield (ref_cov, var_cov)

    @staticmethod
    def bypass_unphased(entries):
        "Filters the given assemblies and annotates the ones without phases so that they bypass coverage calculation"
        unphased_count = 0
        total_var_count = 0
        total_ref_count = 0
        bypassed = []
        for a in entries:
            if a.matches_reference:
                total_ref_count += 1
            else:
                total_var_count += 1
                if not a.phase_ids:
                    unphased_count += 1
                    a.read_coverage = bgexvar.ReadCoverage()
                    a.pair_read_coverage = bgexvar.ReadCoverage()
                    a.bypass_coverage = True
                    bypassed.append(a)
                    if len(bypassed) > 100:
                        bypassed = bypassed[50:]
            yield a
        if unphased_count:
            logging.warning(f"{unphased_count} unphased assemblies encountered out of {total_var_count} variant, {total_ref_count} ref assemblies")
            if unphased_count >= total_var_count:
                for b in bypassed:
                    logging.warning(f"Bypassed phases {b.phase_ids}: {b} ")

    @staticmethod
    def sort_priority_for_align(a):
        """Sort priority for resolving ambiguously mapped reads; higher return value
is higher priority"""
        return (
            # Sort reference first, so we always have reference coverage present
            a.matches_reference,
            # Otherwise, prefer things with less structural variant length
            -abs(a.right_offset - a.left_offset - len(a.seq)),
            # Otherwise, prefer things with longer sequence length
            len(a.seq),
            # Otherwise, prefer things with more supporting reads
            len(a.read_coverage))


    def sort_block_for_align(self, asms):
        """Sort a block for ambiguous read alignment resolution; earlier in the block is higher priority to place reads."""
        asms.sort(key=self.sort_priority_for_align, reverse=True)
        return asms

    @classmethod
    def verify_order(cls, entries):
        """Ensures that incoming entries are in the proper deterministic order, and warns once if they aren't."""
        last = None
        order_errors = 0
        tot_entries = 0
        for e in entries:
            tot_entries += 1
            if last is not None:
                if e < last:
                    if not order_errors:
                        logging.warning(f"Internal error: Saw {last} before {e}")
                    order_errors += 1
                elif not (last < e):
                    if not order_errors:
                        logging.warning(f"Potentially nondeterministic order between {last} and {e}")
                    order_errors += 1
                yield last
            last = e
        if last is not None:
            yield last
        if order_errors:
            logging.warning(f"{order_errors} nondeterministic orderings detected out of {tot_entries} ({order_errors*100./tot_entries:.2f}%)")

    def parse(self, rmap, entries):
        """
        Annotate a set of bgexvar.Assembly entries by editing the VcfEntryInfo
        yields the new bgexvar.Assembly objects
        """
        # Multi-allelic is going to be weird here
        self.readmap = rmap

        if self.phasing:
            entries = self.bypass_unphased(entries)
            entries = bgexvar.join_phases(entries, rmap.max_read_len(), 100)
            if self.max_coverage_alleles:
                entries = LimitAlleles(self.max_coverage_alleles).parse(entries)

        entries = bgexvar.generate_read_cov(rmap, entries,
                                            max_reads_per_entry=self.max_reads_per_entry,
                                            max_coverage_paths=self.max_coverage_paths)

        entries = self.verify_order(entries)

        if self.filter_dup_align:
            entries = bgexvar.filter_dup_align(self.sort_block_for_align, entries)
        if self.ideal_insert:
            entries = bgexvar.place_pair_cov(entries, rm=self.readmap,
                                             min_insert_size=self.min_dist,
                                             max_insert_size=self.max_dist,
                                             ideal_insert_size=self.ideal_insert,
                                             max_ambig=self.placer_max_ambig)
        else:
            entries = bgexvar.generate_pair_cov(rmap, entries, self.min_dist, self.max_dist)


        if self.phasing:
            entries = map(bgexvar.propagate_subassembly_coverage, entries)
            entries = bgexvar.split_phases(entries)

        entries = map(self.with_depths, entries)
        entries = bgexvar.apply_edges(entries, self.annotate_from_edge)
        entries = self.extract_cov_pass_ref(entries)

        for ref, alt in entries:
            if alt:
                yield self.build_vcf(ref, alt)
            else:
                # ref is generated by ADD_REF and should be passed on
                # to the pair coverage generator (within the genotype
                # annotator), so that it can # see pairing data around
                # assemblies.
                yield ref

    @staticmethod
    def get_blank_format():
        """
        Returns the default values
        """
        return AlleleCoverageBase().to_samp_dict()

    # Returns the header information we're creating
    @staticmethod
    def get_header():
        """
        Pull up the header information from AlleleCoverage
        returns list of strings that should be put in the header
        """
        return AlleleCoverageBase.get_header()

    @staticmethod
    def get_format_tags():
        """
        Pull up the tags populated inside an AlleleCoverage
        returns list of tag strings
        """
        return AlleleCoverageBase.get_format_tags()

class LimitAlleles:
    "Processes assemblies to limit the number of simultaneous alleles"

    def __init__(self, max_alleles):
        self.max_alleles = max_alleles
        self.limited_count = 0
        self.sort_count = 0

    @staticmethod
    def sort_priority(a):
        "Returns a larger value for assemblies with a higher priority."
        return (
            # Sort reference first, so we always have reference coverage present
            a.matches_reference,
            # Otherwise, prefer things with more phases present,
            # e.g. more samples which contain this variant.
            len(a.phase_ids),
            # Otherwise, prefer things with more bases, ideally in both reference and sequence.
            (a.right_offset - a.left_offset) + len(a.seq))

    def sort_alleles(self, asms):
        "Sorts assemblies with the highest priority assemblies first"
        #        logging.debug(f"Sorting {len(asms)} alleles from {min(a.left_offset for a in asms)} to {max(a.right_offset for a in asms)}")
        self.sort_count += 1
        asms.sort(key=self.sort_priority, reverse=True)
        return asms

    def on_limited(self, a):
        "Record that an assembly has been allele-limited"
        a.read_coverage = bgexvar.ReadCoverage()
        a.pair_read_coverage = bgexvar.ReadCoverage()
        a.bypass_coverage = True
        a.phase_ids.clear()
        self.limited_count += 1

    def parse(self, entries):
        "Processes the given assemblies"

        entries = bgexvar.limit_alleles(self.max_alleles, self.sort_alleles, self.on_limited, entries)
        total_count = 0
        for a in entries:
            yield a
            total_count += 1
        if self.limited_count:
            logging.info(f"Discarded {self.limited_count} of {total_count} ({self.limited_count * 100 / total_count:.2f}%) assemblies due to allele limits in {self.sort_count} block(s)")
