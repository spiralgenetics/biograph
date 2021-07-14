"""
Genotype variants in a sample
"""
import sys
import json
import argparse
import itertools
import multiprocessing
from io import BytesIO
from collections import namedtuple
import vcf
import biograph as bgsdk
import biograph.internal as bginternal
from biograph.internal.kPCMP import edit_vcf_header
from biograph.utils import genotyper, ConsumerPool
import biograph.tools.log as log
import progressbar

class AnnoEntry:
    """
    Edits a vcf entry in-place
    """

    def __init__(self, entry, sample=None, gt_replace=False):
        """
        Hold an entry for annotating
        Sample will only place PGT and such in given sample
        Replace GT with PGT when gt_replace
        """
        self.entry = entry
        self.sample = sample
        self.gt_replace = gt_replace
        self.__new_cand = []
        self.__alleles = []
        self.__update_data()

    def __update_data(self):
        """
        Places the new annotations in the entry
        """
        if self.entry.FORMAT is None:
            self.entry.FORMAT = "PR"
        fmt_data = self.entry.FORMAT.split(':')
        if "PR" not in fmt_data:
            self.entry.add_format("PR")
        if "PA" not in fmt_data:
            self.entry.add_format("PA")
        if "PGT" not in fmt_data:
            self.entry.add_format("PGT")
        if "PGQ" not in fmt_data:
            self.entry.add_format("PGQ")

    def add_allele(self, genoid, paired_cov, unpaired_cov):
        """
        Collect the allele coverage information
        """
        self.__alleles.append((genoid, paired_cov, unpaired_cov))

    def iter_alleles(self):
        """
        Iterates through all pairs of alleles
        """
        for allele1, allele2 in itertools.combinations(self.__alleles, 2):
            yield allele1, allele2

    def add_data(self, new_pr, new_pa, new_pgt, new_pgq):
        """
        a potential genotype canidate
        """
        self.__new_cand.append([new_pr, new_pa, new_pgt, new_pgq])

    def update_entry(self):
        """
        edit the entry with the new genotype candidate
        """
        if not self.__new_cand:
            n_data = ('.', '.', '.', '.')
        else:
            log.debug("all candidates for %s", str(self.entry))
            log.debug("new_pr, new_pa, new_pgt, new_pgq")
            for i in self.__new_cand:
                log.debug(i)
            self.__new_cand.sort(cmp=lambda x, y: int(round(x[3] - y[3])), reverse=True)
            n_data = tuple(self.__new_cand[0])

        if self.sample is None:
            v_sample = self.entry.samples[0]
        else:
            v_sample = self.entry.genotype(self.sample)

        n_call_data = vcf.model.make_calldata_tuple(v_sample.data._fields + ("PR", "PA", "PGT", "PGQ"))

        v_sample.data += n_data
        v_sample.data = n_call_data(*v_sample.data)

    def add_genotype(self, vgraph):
        """
        Creates all the annotations around genotype for the vcf entry
        """
        #for each allele in the entry, I need to use vgraph to get the coverages
        ref_cov = vgraph.get_allele_coverage(self.entry.start, self.entry.end)
        self.add_allele("0", *ref_cov)

        for genoid, allele in enumerate(self.entry.ALT):
            paired, unpaired = vgraph.get_allele_coverage(self.entry.start, self.entry.end, allele.sequence)
            self.add_allele(str(genoid + 1), paired, unpaired)

        for allele1, allele2 in self.iter_alleles():
            #coverage over the two alleles
            a_cov = [allele1[1], allele2[1]]
            log.debug("checking alleles %s %s", str(allele1), str(allele2))
            # Unoptimal to hardcode these
            #if sum(a_cov) < 5:  #Need to flag as used unpaired..
                #log.info("not enough cov %s\n%s", str(a_cov), str(self.entry))
                #a_cov[0] += int(round(.5 * allele1[2]))
                #a_cov[1] += int(round(.5 * allele2[2]))
            # imagine C, ALT=[CATTTT, T], that's 2/2
            # allele #1 is going to have zero ref or alt coverage.
            # and therefore must be ignored - throwing this
            # into bgsdk.genotyper gives it a 'GQ' of 20... horrible
            new_gt, new_gq = genotyper(
                a_cov[0] + a_cov[1],
                a_cov[1],
                priors=[["{0}/{0}".format(allele1[0]), 0.01], [allele1[0] + '/' + allele2[0], 0.50],
                        ['{0}/{0}'.format(allele2[0]), .99]])
            self.add_data(a_cov[0], a_cov[1], new_gt, new_gq)

        self.update_entry()

class VarGraph:
    """ Makes it easy to look up nodes/edges of vcf records in a vgraph """

    def __init__(self, bg, ref, min_insert=200, max_insert=1000):
        """
        Given a biograph, a reference range, and a vcf,
        build a VarGraph and wrap it's functionality
        """
        # Need to hide all these probably
        self.bg = bg
        self.ref = ref
        self.min_insert = min_insert
        self.max_insert = max_insert

        self.vgraph = bginternal.VarGraph(self.ref.sequence, min_insert, max_insert)
        self.__ref_alleles_lookup = {}
        self.__alt_alleles_lookup = {}
        self.var_count = 0

    def __shift(self, entry):
        """
        Without trimming the record,
        move the vcf variant  to it's relative position within the trace span
        this is just a helper for load_vcf
        """
        for alt_allele in entry.ALT:
            if not hasattr(alt_allele, "sequence"):
                continue
            start = entry.start - self.ref.start
            end = entry.end - self.ref.start
            alt_seq = bgsdk.Sequence(alt_allele.sequence.replace("N", "A"))
            yield start, end, alt_seq

    def __shift_trim(self, entry):
        """
        Move a vcf_record to it's relative position on the
        reference_range
        """
        for alt_allele in entry.ALT:
            if not hasattr(alt_allele, "sequence"):
                continue
            start = entry.start - self.ref.start
            end = entry.end - self.ref.start
            alt_seq = alt_allele.sequence
            if "N" in alt_seq:
                log.warning("Entry %s has N in ALT allele. Turning to 'A'", str(entry))
                alt_seq = alt_seq.replace('N', 'A')
            if len(alt_seq) == len(entry.REF):
                alt_seq = bgsdk.Sequence(alt_seq)
            elif len(entry.ALT[0]) == 1:  # presumed deletion
                if str(entry.ALT[0])[0] != entry.REF[0]:  # except the repl types
                    alt_seq = bgsdk.Sequence(alt_seq)
                else:
                    # Get rid of anchor base, no alternate sequence
                    start += 1
                    alt_seq = bgsdk.Sequence("")
            else:  # presumed insertion or repl
                # Get rid of anchor base in alt seq as well
                start += 1
                alt_seq = bgsdk.Sequence(alt_seq[1:])
            yield start, end, alt_seq

    def load_vcf(self, vcf_file):
        """
        Add the vcf entries over the vargraph region into the graph
        """
        for entry in vcf_file.fetch(self.ref.chromosome, self.ref.start, self.ref.end):
            if entry.start <= self.ref.start or entry.end >= self.ref.end:
                #We have a gap spanning SV. skip it
                continue
            for start, end, alt_seq in self.__shift(entry):
                if alt_seq is not None:  # alt_seq is None for 'weird' variants
                    try:
                        self.add_variant(start, end, alt_seq)
                    except RuntimeError:
                        log.critical("Bad variant coords (%d, %d) in (%d, %d) @ %s", start, end, self.ref.start,
                                     self.ref.end, str(entry))

        # if you don't load any variants, you don't need to trace it...
        log.debug("%d variants in region %s:%d-%d", self.var_count, self.ref.chromosome, self.ref.start,
                  self.ref.end)

    def add_variant(self, start, end, alt_seq):
        """
        Wrapper around vgraph
        """
        self.var_count += 1
        self.vgraph.add_variant(start, end, alt_seq)

    def __build_lookup(self):
        """
        Make node lookups for alternate and reference alleles.
        """
        for node in self.vgraph.get_nodes():
            if node.is_ref:
                #Edge that goes into this node from reference
                for edge in node.upstream:
                    if edge.up_node.is_ref and edge.dn_node.is_ref and \
                      edge.up_node.end == node.start:
                        self.__ref_alleles_lookup[node.start] = edge
                for edge in node.downstream:
                    if edge.up_node.is_ref and edge.dn_node.is_ref and \
                      edge.dn_node.start == node.end:
                        self.__ref_alleles_lookup[node.end] = edge
            else:
                self.__alt_alleles_lookup[(node.start, node.end, node.seq)] = (node.upstream[0], node.downstream[0])

    def trace(self, sample=None):
        """
        simply calls trace
        """
        if self.var_count == 0:
            return
        if sample is not None:
            self.bg.set_readmap(sample)
        self.vgraph.trace(self.bg)
        self.__build_lookup()

    def get_allele_coverage(self, start, end, sequence=None):
        """
        Gets the coverage of an allele in the variant graph.
        If the sequence is None, return a reference allele
        Returns (pairedcov, unpairedcov)
        """
        start -= self.ref.start
        end -= self.ref.start
        if sequence is None:
            upedge = self.__ref_alleles_lookup[start]
            dnedge = self.__ref_alleles_lookup[end]
        else:
            if isinstance(sequence, str):
                sequence = bgsdk.Sequence(sequence.replace("N", "A"))
            upedge, dnedge = self.__alt_alleles_lookup[(start, end, sequence)]
        paired = int(round((upedge.paired + dnedge.paired))/2)
        unpaired = int(round((upedge.unpaired + dnedge.unpaired))/2)
        return paired, unpaired

    def dump_graph(self):
        """
        Return String representation of Graph
        """
        return self.vgraph.dump_graph()

class PcmpTask:
    """
    This is one run of the program - I can thread these
    """

    def __init__(self, bg_file, ref_file, var_file, region=None, \
                 sample=None, gt_replace=False, min_insert=200, max_insert=1000):
        """
        Main runner
        On a BioGraph file
        and a Reference
        genotype the variants in var_file
        """
        self.bg_file = bg_file
        self.ref_file = ref_file
        self.var_file = var_file
        self.region = region
        self.sample = sample
        self.gt_replace = gt_replace
        self.min_insert = min_insert
        self.max_insert = max_insert
        # This needs to be changed to a vcf StringIO or something
        # And we need a vcf header editor / sample data changer thing
        self.annotated_vars = None
        self.name = "%s:%d-%d" % tuple(self.region)

    def collect_alleles(self, vgraph, ref_region):
        """
        Given a traced variant graph and its reference region,
        create an output vcf.
        For every variant, goin var_lookup, get the pieces it needs to be genotyped
        Change this to write genotypes - or something.
        Returns a full vcf in StringIO
        annotate on a single sample with sample
        if gt_replace: When GT is blank, replace with PGT
        """
        vcf_file = vcf.Reader(filename=self.var_file)
        edit_vcf_header(vcf_file)

        self.annotated_vars = BytesIO()
        ret_vcf = vcf.Writer(self.annotated_vars, vcf_file)
        try:
            for entry in vcf_file.fetch(ref_region.chromosome, ref_region.start, ref_region.end):
                if not (entry.start >= ref_region.start and entry.start <= ref_region.end):
                    continue  #prevent redundant entries.
                anno = AnnoEntry(entry, self.sample, self.gt_replace)
                if entry.start <= ref_region.start or entry.end >= ref_region.end: #already handled?
                    log.debug("Entry %s spans a contig. Cannot be properly traced.", str(entry))
                    anno.update_entry()
                else:
                    anno.add_genotype(vgraph)
                ret_vcf.write_record(anno.entry)
        except ValueError:
            pass  #Tabix problems sometimes
        ret_vcf.flush()

    def __call__(self):
        """
        Runs all the steps to perform pcmp.
        self.all_variants is updated in place
        """
        m_bg = bgsdk.BioGraph(self.bg_file, bgsdk.CacheStrategy.MMAP)

        m_ref = bgsdk.Reference(self.ref_file)
        try:
            ref_region = m_ref.make_range(self.region[0], self.region[1], self.region[2])
        except RuntimeError:
            log.error("Reference range %s has Ns, cannot parse", self.name)
            return

        # Easier to handle
        log.debug("Running pcmp on %s", self.name)
        vlook = VarGraph(m_bg, ref_region, self.min_insert, self.max_insert)
        m_vars = vcf.Reader(filename=self.var_file)
        log.debug("Loading variants on %s", self.name)
        vlook.load_vcf(m_vars)
        log.debug("Beginning trace on %s", self.name)
        vlook.trace(self.sample)
        # Make a new vcf with annos
        # For each variant, var_graph.whatever
        log.debug("End trace on %s", self.name)
        log.debug("Begin collection on %s", self.name)
        self.collect_alleles(vlook, ref_region)
        log.debug("%s dump\n %s", self.name, vlook.dump_graph())
        log.debug("Finished %s", self.name)


class Pcmp:  # pylint: disable=too-many-instance-attributes
    """
    Wrapper around the entire pcmp program
    """

    def __init__(self, # pylint: disable=too-many-arguments
                 bg_fn,
                 ref_fn,
                 vcf_fn,
                 threads=1,
                 out="out.vcf",
                 min_insert=200,
                 max_insert=1000,
                 sample=None,
                 region=None,
                 gt_replace=False):
        """
        Hold the variables
        """
        self.bg_fn = bg_fn
        self.ref_fn = ref_fn
        self.vcf_fn = vcf_fn
        self.m_ref = bgsdk.Reference(ref_fn)
        self.vcf_file = vcf.Reader(filename=vcf_fn)
        self.out_file = out
        self.threads = threads
        self.min_insert = min_insert
        self.max_insert = max_insert
        self.sample = sample
        self.gt_replace = gt_replace
        self.regions = self.__parse_regions(region)

    def __parse_regions(self, region):
        """
        Create the individual regions over which this run is going to parse
        """
        FakeSC = namedtuple("fake_supercontig", ["chromosome", "start", "end"])
        ret = []
        if region is None:
            return self.m_ref.supercontigs

        try:
            chrom, pos = region.split(':')
            start, end = [int(x) for x in pos.split("-")]
            for i in self.m_ref.find_ranges(chrom, start, end):
                n_chrom = i.chromosome
                n_start = i.start
                n_end = i.end
                #Problem when --region is within an N
                if n_start < start:
                    n_start = start
                if n_end > end:
                    n_end = end
                ret.append(FakeSC(n_chrom, n_start, n_end))
        except ValueError:
            log.error("Problem parsing region %s - expected chrom:start-end", region)
            exit(1)

        return ret

    def make_range_regions(self):
        """
        Parse all the entires of a vcf file and subset the regions
        returns a generator
        """
        var_cnt = 0
        sub_reg_cnt = 0
        contig_cnt = 0
        for ref_range in self.regions:
            contig_cnt += 1
            # Wrap the region
            cur_start = None
            cur_end = None
            try:
                for entry in self.vcf_file.fetch(ref_range.chromosome, ref_range.start, ref_range.end):
                    var_cnt += 1
                    if cur_start is None:
                        cur_start = max(ref_range.start, entry.start - self.max_insert)
                        cur_end = min(entry.end + self.max_insert, ref_range.end)
                    #elif entry.start <= cur_end: #only insert. not double insert.
                    elif (entry.start - self.max_insert) <= cur_end:
                        # Expand this region if there's overlap
                        cur_end = min(entry.end + self.max_insert, ref_range.end)
                    else:
                        # Make a new region, start recording the nextone
                        yield (ref_range.chromosome, cur_start, cur_end)
                        sub_reg_cnt += 1
                        cur_start = max(ref_range.start, entry.start - self.max_insert)
                        cur_end = min(entry.end + self.max_insert, ref_range.end)
            except ValueError:
                pass #vcf fetching problems sometimes
            # Pick up the last region
            if cur_start is not None:
                yield (ref_range.chromosome, cur_start, cur_end)
            sub_reg_cnt += 1
        log.info("Loaded %d variants in %d trace regions over %d contigs", var_cnt, sub_reg_cnt, contig_cnt)

    def do_pcmp(self):
        """
        The full runner
        """
        procs = ConsumerPool(self.threads)
        procs.start_pool()
        num_regions = 0
        log.info("Creating trace regions")
        try:
            for sub_region in self.make_range_regions():
                procs.put_task(
                    PcmpTask(
                        bg_file=self.bg_fn,
                        ref_file=self.ref_fn,
                        var_file=self.vcf_fn,
                        region=sub_region,
                        sample=self.sample,
                        gt_replace=self.gt_replace,
                        min_insert=self.min_insert,
                        max_insert=self.max_insert)
                )
                num_regions += 1
        #many possiblities..
        except IOError as e:
            log.error("Problem running pcmp. %s", str(e))
        log.info("%d regions loaded for parsing", num_regions)
        procs.put_poison()
        tmp_vcf = vcf.Reader(filename=self.vcf_fn)
        edit_vcf_header(tmp_vcf)
        var_out = vcf.Writer(open(self.out_file, 'w'), tmp_vcf)
        pgbar = progressbar.ProgressBar(redirect_stdout=True, max_value=num_regions, widgets=[
            ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(num_regions), '] ',
            progressbar.Bar(), ' (', progressbar.ETA(), ') ',
        ])
        cnt = 0
        pgbar.update(cnt)
        #KeyboardInterrupt catch here to kill procs
        for pcmp_result in procs.get_tasks():
            cnt += 1
            pgbar.update(cnt)
            if pcmp_result.annotated_vars is not None:
                pcmp_result.annotated_vars.seek(0)
                tmp_vcf = vcf.Reader(pcmp_result.annotated_vars)
                for entry in tmp_vcf:
                    var_out.write_record(entry)
            else:
                log.error("Error in variants over region %s", str(pcmp_result.region))
        var_out.close()


def parse_args(args):
    """ Make pretty arguments """
    parser = argparse.ArgumentParser(description='Genotype variants in a VCF')
    parser.add_argument(
        "-b", "--biograph", metavar="BG", required=True, help="Merged BioGraph file containing individuals")
    parser.add_argument(
        "-v", "--variants", metavar="VCF", required=True, help="The Anchored Assembly vcf containing all samples")
    parser.add_argument("-r", "--reference", metavar="REF", required=True, help="Reference genome folder")
    parser.add_argument("-R", "--region", default=None, help="Region of the reference to pcmp chr:start-end")
    parser.add_argument("-s", "--sample", default=None, help="Sample in merged BioGraph to use(%(default)s)")
    parser.add_argument("-o", "--output", metavar="OUT", default="output.vcf", help="Annotated vcf file output.")
    parser.add_argument(
        "-g", "--gt-replace", default=False, help="Replace existing GT fields with the PGT (%(default)s)")
    parser.add_argument("-m", "--min-insert", default=200, type=int, help="Minimum insert size to consider paired (%(default)s)")
    parser.add_argument("-M", "--max-insert", default=1000, type=int, help="Maximum insert size to consider paired (%(default)s)")
    parser.add_argument(
        "-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true", help="Verbose logging")
    args = parser.parse_args(args)
    log.setup_logging(args.debug)
    log.debug("Params:\n%s", json.dumps(vars(args), indent=4))

    return args


def vpcmp_run(args):
    """ Setup and run pcmp """
    args = parse_args(args)
    log.info("Pre-caching BioGraph")  #do I need to do this? probably.
    tmp_bg = bgsdk.BioGraph(args.biograph, bgsdk.CacheStrategy.MMAPCACHE) # pylint: disable=unused-variable
    #Could infer the insert size here and set those params
    my_p = Pcmp(
        bg_fn=args.biograph,
        ref_fn=args.reference,
        vcf_fn=args.variants,
        threads=args.threads,
        out=args.output,
        min_insert=args.min_insert,
        max_insert=args.max_insert,
        region=args.region,
        sample=args.sample,
        gt_replace=args.gt_replace)
    my_p.do_pcmp()
    log.info("Finished")


if __name__ == '__main__':
    vpcmp_run(sys.argv[1:])
