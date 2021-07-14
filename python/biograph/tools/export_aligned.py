"""
Calculates coverage over an input VCF and Outputs aligned reads in BAM format.
"""
import sys
import json
import argparse
import multiprocessing

from biograph.tools import ParallelRegions
import biograph.coverage as bganno
import biograph.variants as bgexvar
import biograph.tools.log as log

class AlignedWriter(multiprocessing.Process):
    """Combines SAM lines from subprocesses to an output file"""
    def __init__(self, output_name):
        super().__init__()
        self.output_name = output_name
        self.queue = multiprocessing.Queue()
        self.outf = None
        self.sam_output = None

    def run(self):
        self.outf = open(self.output_name, "wb")
        while True:
            data = self.queue.get()
            if data is None:
                print(f"Done Data")
                break
            self.outf.write(data)
        print(f"Closing")
        self.outf.close()
        print(f"Done closing")

class ExportAligned(ParallelRegions):
    """Export aligned reads in BAM format"""
    def __init__(self):
        super().__init__("export_aligned")
        self.output_bam = None
        self.output_queue = None
        self.args = None

    def close_output(self):
        "Closes output SAM"
        self.output_queue.put(None)
        self.output_queue.close()

    def parse_args(self, clargs):
        """Parses command line arguments"""
        parser = argparse.ArgumentParser(prog="export_aligned", description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-o", "--output", metavar="OUT", default="output.bam",
                            help="Resultant aligned reads.")
        parser.add_argument("--skipref-on-insert", default=False, action="store_true",
                            help="Add a '1N' cigar operation before reads starting inside insertions.  This helps some alignment viewers display the insertions properly.")
        parser.add_argument("--emit-all-mates", default=False, action="store_true",
                            help="Make sure mates are emitted for all paired reads, even if they're not mapped to any variants present")
        parser.add_argument("--unaligned", default=False, action="store_true",
                            help="Don't bother aligning; just output all reads as unmapped reads")
        parser.add_argument("--debug", action="store_true", help="Verbose logging")
        parser.add_argument("--ideal-insert", type=int, default=None,
                            help="Only export pairs, and only export each pair once, at the closest insert size to ideal.")
        parser.add_argument("--max-insert", type=int, default=None,
                            help="Minimum insert size")
        parser.add_argument("--min-insert", type=int, default=None,
                            help="Maximum insert size")
        self.setup_funcs.append(self.setup_export)
        self.add_reference_arguments(parser)
        self.add_vcf_arguments(parser)
        self.add_biograph_arguments(parser)
        self.add_readmap_arguments(parser)
        self.add_read_coverage_arguments(parser)
        return parser.parse_args(clargs)

    def setup_export(self, args):
        """Opens files for export"""
        log.setup_logging(args.debug)
        log.debug("Params:\n%s", json.dumps(vars(args), indent=4))
        self.output_bam = args.output
        self.args = args

        if args.ideal_insert or args.max_insert or args.min_insert:
            if not (args.max_insert and args.min_insert):
                raise RuntimeError("Must specify both maximum and minimum insert size if using pairing data")

    def open_output(self, queue):
        "Starts writing output SAM through the queue to the AlignedWriter"
        self.output_queue = queue
        sam_output = bganno.SamOutput(self.readmap)
        sam_output.add_reference(self.reference)
        for line in sam_output.get_header():
            self.output_queue.put(line)

    @staticmethod
    def sort_priority(a):
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

    def sort_block(self, asms):
        """Sort a block for ambiguous read alignment resolution; earlier in the block is higher priority to place reads."""
        asms.sort(key=self.sort_priority, reverse=True)
        return asms

    @staticmethod
    def copy_pair_to_read(entries):
        """Copies from pair_read_coverage to read_coverage after calcultating pairing data"""
        for a in entries:
            a.read_coverage = a.pair_read_coverage
            yield a

    def process_region(self, region):
        """Exports aligned reads in the given region"""
        entries = self.fetch_vcf(region)
        entries = self.vcf_to_assemblies(entries, save_vcf_line=False)
        entries = self.add_reference(region, entries)
        entries = self.join_phases(entries)
        entries = self.generate_read_coverage(entries)
        entries = self.split_phases(entries)
        if self.args.ideal_insert:
            entries = bgexvar.place_pair_cov(entries, rm=self.readmap,
                                             ideal_insert_size=self.args.ideal_insert)
            entries = self.copy_pair_to_read(entries)
        elif self.args.max_insert and self.args.min_insert:
            entries = bgexvar.generate_pair_cov(input=entries, rm=self.readmap,
                                                min_insert_size=self.args.min_insert, max_insert_size=self.args.max_insert)
            entries = self.copy_pair_to_read(entries)
        elif not self.args.unaligned:
            entries = bgexvar.filter_dup_align(self.sort_block, entries)

        entries = self.output_sam(region, entries)
        for _ in entries:
            pass

    def output_sam(self, region, entries):
        """Exports aligned reads from the given assemblies"""
        ctg, _, _ = region
        sam_output = bganno.SamOutput(self.readmap, emit_all_mates=self.args.emit_all_mates)
        if self.args.unaligned:
            yield from sam_output.output_unaligned(self.output_queue, ctg, entries)
        else:
            sam_output.add_reference(self.reference)
            yield from sam_output.output_aligned(self.output_queue, ctg, entries)

def main(clargs):
    """
    Main method to run export aligned reads
    """


    exporter = ExportAligned()
    args = exporter.parse_args(clargs)
    writer = AlignedWriter(args.output)
    writer.start()
    exporter.setup(args)
    exporter.open_output(writer.queue)
    exporter.process_regions()
    log.info("Finished processing regions")
    exporter.close_output()
    writer.join()

if __name__ == '__main__':
    main(sys.argv[1:])
