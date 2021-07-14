"""Provides spuport for processing  a set of reference regions in parallel.

Provides support for command line arguments, reading VCFs and
biographs, converting VCF lines to assemblies, and calculating
coverage.

"""

import multiprocessing
from types import SimpleNamespace
import logging
import sys
import os
from setproctitle import setproctitle  # pylint: disable=no-name-in-module

import biograph.variants as bgexvar
import biograph.coverage as bganno
import biograph

import tabix

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

class ParallelRegions:
    """Represents a processor to process regions of reference in parallel.

Users can subclass this class and customize "process_region"

"""

    instance = None

    def __init__(self, progname):
        self.reference = None
        self.biograph = None
        self.readmap = None
        self.threads = None
        self.regions = None
        self.variant_file = None
        self.coverage = None
        self.setup_funcs = []
        self.progname = progname

    def add_reference_arguments(self, parser):
        """Add basic relating to references, regions, and threads"""
        parser.add_argument("-r", "--reference", metavar="REF", required=True,
                            help="Reference genome folder")
        parser.add_argument("-R", "--regions", "--bed", default=None,
                            help="Bed file of regions to process")
        parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                            help="Number of threads to use (%(default)s)")
        self.setup_funcs.append(self.setup_reference)

    def setup_reference(self, args):
        """Configures the reference after arguments have been parsed."""
        self.reference = biograph.Reference(args.reference)
        self.threads = args.threads

        self.regions = []
        if args.regions is not None:
            with open(args.regions, 'r') as fh:
                for line in fh:
                    data = line.strip().split('\t')
                    self.regions.append((data[0], int(data[1]), int(data[2])))
        else:
            for ctg in self.reference.scaffolds:
                self.regions.append((ctg, 0, int(self.reference.scaffold_lens[ctg])))

    def add_vcf_arguments(self, parser):
        """Add arguments relating to reading VCF files"""
        parser.add_argument("-v", "--variants", metavar="VCF", required=True,
                            help="Input VCF")
        self.setup_funcs.append(self.setup_vcf)

    def setup_vcf(self, args):
        """Configures VCF file reading after arguments have been parsed."""
        if not args.variants.endswith(".vcf.gz") or not os.path.exists(args.variants + ".tbi"):
            logging.error("Variants file must be compressed and indexed")
            sys.exit(1)
        self.variant_file = args.variants

    def add_biograph_arguments(self, parser):
        """Add arguments relating to reading BioGraph files"""
        parser.add_argument("-b", "--biograph", metavar="BG", required=True,
                            help="Merged BioGraph file containing individuals")
        parser.add_argument("--cache", default=False, action="store_true",
                            help="Attempt to cache as much as possible in RAM")
        self.setup_funcs.append(self.setup_biograph)

    def setup_biograph(self, args):
        """Opens BioGraph file after arguments hvae been parsed."""
        self.biograph = biograph.BioGraph(args.biograph,
                                          biograph.CacheStrategy.RAM if args.cache
                                          else biograph.CacheStrategy.MMAPCACHE)

    def add_readmap_arguments(self, parser):
        """Configures arguments relating to using a readmap"""
        parser.add_argument("-s", "--sample", default=None,
                            help="Sample in merged BioGraph to use (%(default)s)")
        self.setup_funcs.append(self.setup_readmap)

    def setup_readmap(self, args):
        """Opens a readmap after arguments have been parsed"""
        if not self.biograph:
            raise RuntimeError("Readmaps require a Biograph")
        samples = self.biograph.metadata.samples
        if args.sample is None:
            if len(samples) == 1:
                args.sample = next(iter(samples))
        elif args.sample not in samples:
            raise KeyError("Sample %s not present in BioGraph" % args.sample)
        self.readmap = self.biograph.open_readmap(args.sample)

    def add_read_coverage_arguments(self, parser):
        """Adds arguments relating to read coverage calculation"""
        parser.add_argument("--max-read-cov-paths", type=int, default=3000,
                            help="If nonzero, limit number of parallel coverage paths to trace through the seqset when calculating coverage (%(default)s)")
        parser.add_argument("--phasing", type=bool, default=False, help="Skip paths counterindicated by phasing information")
        parser.add_argument("--max-paths", type=int, default=100,
                            help="If nonzero, limit number of paths through VCF entries to the given value.")
        parser.add_argument("--max-reads-per-entry", type=int, default=0,
                            help="If nonzero, maximum number of reads per seqset entry to use; additional reads are ignored (%(default)s)")
        self.setup_funcs.append(self.setup_read_coverage)

    def setup_read_coverage(self, args):
        """Saves arguments relating to calculating read coverage"""
        if not self.readmap:
            raise RuntimeError("Coverage requires a readmap")
        self.coverage = SimpleNamespace(
            max_read_cov_paths=args.max_read_cov_paths,
            phasing=args.phasing,
            max_paths=args.max_paths,
            max_reads_per_entry=args.max_reads_per_entry)

    def setup(self, args):
        """Configures processing after arguments have been parsed"""
        for setup_func in self.setup_funcs:
            setup_func(args)

    def fetch_vcf(self, region):
        """Generator which returns all vcf lines in the given region"""
        ctg, start, stop = region
        tb = tabix.open(self.variant_file)
        try:
            for line in tb.query(ctg, start, stop):
                yield line
        except tabix.TabixError as e:
            logging.error(f"Tabix encountered error: {e}")

    def vcf_to_assemblies(self, entries, *, save_vcf_line=True):
        """Generator which converts the given vcf lines into BioGraph assemblies"""
        assembly_id = 0
        did_warn = False
        for vcf_line in entries:
            phase_ids = None
            if self.coverage and self.coverage.phasing:
                formats = vcf_line[8].split(':')
                phasing_format_offset = formats.index("PI")
                phase_ids = bgexvar.PhaseSet.from_format_fields(phasing_format_offset, vcf_line[9:])
            start = int(vcf_line[1]) - 1
            stop = start + len(vcf_line[3])
            for alt in vcf_line[4].split(","):
                assembly_id += 1
                try:
                    a = bgexvar.Assembly(start, stop, alt, assembly_id)
                except TypeError:
                    if not did_warn:
                        did_warn = True
                        logging.warning(f"Couldn't create assembly from {start} to {stop} with alt={alt}")
                    continue
                if save_vcf_line:
                    a.vcf_line = vcf_line
                if phase_ids:
                    a.phase_ids = phase_ids
                yield a

    def add_reference(self, region, entries):
        """Generator which adds reference assemblies to the given assembly stream"""
        ctg, _, _ = region
        entries = bgexvar.trim_ref(self.reference, ctg, entries)
        entries = bgexvar.add_ref_assemblies(self.reference, ctg, entries, self.readmap.max_read_len())
        return entries

    def join_phases(self, entries):
        """Generator which joins like phases within the given assembly stream if phasing is enabled"""
        if self.coverage.phasing:
            entries = bganno.PhaseConflictResolver().parse(None, entries)
            entries = bypass_unphased(entries)
            entries = bgexvar.join_phases(entries, self.readmap.max_read_len(), 100000)
        return entries

    def generate_read_coverage(self, entries):
        """Generator which calculates coverage for the given assemblies"""
        if self.coverage.max_paths:
            entries = LimitAlleles(self.coverage.max_paths).parse(entries)
        entries = bgexvar.generate_read_cov(self.readmap,
                                            entries,
                                            max_reads_per_entry=self.coverage.max_reads_per_entry,
                                            max_coverage_paths=self.coverage.max_read_cov_paths)
        return entries

    def split_phases(self, entries):
        """Generator which splits phased assemblies"""
        if self.coverage.phasing:
            entries = map(bgexvar.propagate_subassembly_coverage, entries)
            entries = bgexvar.split_phases(entries)
        return entries

    def log_and_process_region(self, region):
        """Processes the given region, logging start and stop times"""
        ctg, start, stop = region
        logging.info(f"Starting {ctg}:{start}-{stop}")
        setproctitle(f"{self.progname} worker: {ctg}:{start}-{stop}")
        output = self.process_region(region)
        logging.info(f"Done {ctg}:{start}-{stop}")
        return output

    @staticmethod
    def process_region(region):
        """Processes the given region.  Subclasses should override this."""
        raise RuntimeError("Subclass should override process_region")

    @classmethod
    def process_region_in_instance(cls, region):
        """Helper function so that we don't have to pickle a ParallelRegion"""
        return ParallelRegions.instance.log_and_process_region(region)

    def process_regions(self):
        """Processes all configured regions in parallel"""
        ParallelRegions.instance = self
        if self.threads > 1:
            with multiprocessing.Pool(self.threads, maxtasksperchild=1) as pool:
                results = list(pool.imap_unordered(
                    ParallelRegions.process_region_in_instance,
                    self.regions))
                pool.close()
                pool.join()
        else:
            results = list(map(ParallelRegions.process_region_in_instance, self.regions))
        ParallelRegions.instance = None
        return results
