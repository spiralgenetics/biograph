"""
Import formats: Test biograph importing input data in various formats
"""

from __future__ import print_function

import unittest
import os
from subprocess import call

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)
from python.functest.utils.defaults import (
    GOLDEN_DIR
)
import biograph

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class ImportFormatTestCases(unittest.TestCase):
    """Test formats"""

    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit code matches (zero by default) """
        self.assertEqual(code, call(cmd, shell=True), "{0} did not exit {1}".format(cmd, code))

    def index_filenames(self, base_filename):
        """Returns a list of all index files for the given data file that exist"""

        index_filenames = []
        data_filenames = []
        for data_suffix in ["", ".bam", ".cram", ".sam"]:
            if data_suffix:
                data_filename = base_filename + data_suffix
                if os.path.exists(data_filename):
                    data_filenames.append(data_filename)
            for index_suffix in [".bai", ".crai", ".csi", ".tbi"]:
                index_filename = base_filename + data_suffix + index_suffix
                if os.path.exists(index_filename):
                    index_filenames.append(index_filename)
        self.assertEqual(len(data_filenames), 1, msg=data_filenames)
        return index_filenames

    def import_file(self, out_name, *args, pipe_from=None, data_format=None):
        """
self.import_file(name, filespec1, filespec2, ...) Imports one or more files
and returns the resulting pairing statistics.

Each filespec may either be a filename or a pair of (read data, pair data)
        """

        self.assertTrue(args)

        # Reference is important, since it is used to decode cram.
        reference = "datasets/reference/e_coli_k12_ASM584v1/"

        cmd = "bgbinary create"
        for filespec in args:
            cmd += " "
            if isinstance(filespec, tuple):
                read_filename, pair_filename = filespec
                cmd += f" --in {read_filename} --pair {pair_filename}"
            else:
                read_filename = filespec
                cmd += f" --in {read_filename}"

        test_biograph = f"{self.data_dir}/{out_name}.bg"
        cmd += f" --ref {reference} --out {test_biograph}"

        # Disable read correction so we get exactly the same file contents.
        cmd += " --min-kmer-count 1"

        # These are all small datasets; don't use much RAM
        cmd += " --max-mem 4"

        # biograph is unhappy if it doesn't get enough RAM
        cmd += f" --tmp {self.data_dir}/tmp"

        if data_format:
            cmd += f" --format {data_format}"
        if pipe_from:
            self.assertTrue(data_format)
            cmd = f"cat {pipe_from} | " + cmd
        self.check(cmd)

        bg = biograph.BioGraph(test_biograph, biograph.CacheStrategy.MMAP)
        rm = bg.open_readmap()
        stats = rm.get_pair_stats()

        print(f"Loaded biograph {out_name}.bg:")
        print(f"  Seqset entries: {bg.seqset.size()}")
        rm = bg.open_readmap()
        print(f"  Readmap entries: {rm.size()}")
        ps = rm.get_pair_stats()
        print(f"  Paired reads: {ps.paired_reads}")
        print(f"  Unpaired reads: {ps.unpaired_reads}")

        return stats

    def test_unpaired_fastq(self):
        """Test unpaired fastq"""

        stats = self.import_file("unpaired_fastq",
                                 "datasets/bams/pairing_test/pair_test1.fq")
        self.assertEqual(stats.paired_reads, 0)
        self.assertEqual(stats.unpaired_reads, 25050)

    def test_stdin_fastq(self):
        """Test fastq supplied on standard input"""

        stats = self.import_file("stdin_fastq",
                                 "-",
                                 pipe_from="datasets/bams/pairing_test/pair_test1.fq",
                                 data_format="fastq")
        self.assertEqual(stats.paired_reads, 0)
        self.assertEqual(stats.unpaired_reads, 25050)

    def test_paired_fastq(self):
        """Test unpaired fastq"""

        orig_fq1 = "datasets/bams/pairing_test/pair_test1.fq"
        orig_fq2 = "datasets/bams/pairing_test/pair_test2.fq"
        truncated_fq1 = f"{self.data_dir}/pair_test1-truncated.fq"
        self.check(f"ls -l {self.data_dir}")
        self.check(f"head -n 50000 {orig_fq1} > {truncated_fq1}")

        stats = self.import_file("paired_fastq", (truncated_fq1, orig_fq2))
        self.assertEqual(stats.paired_reads, 25000)
        self.assertEqual(stats.unpaired_reads, 12550)

    def test_unindexed_bam(self):
        """Test unindexed BAM file"""

        base_filename = f"{GOLDEN_DIR}/seqset/hiv_test"
        self.assertEqual(self.index_filenames(base_filename), [],
                         msg="We're trying to test without indexes here")
        bam_filename = base_filename + ".bam"

        stats = self.import_file("unindexed_bam", bam_filename)
        self.assertEqual(stats.paired_reads, 998)
        self.assertEqual(stats.unpaired_reads, 1)

    def test_stdin_bam(self):
        """Test BAM file supplied on standard input"""

        stats = self.import_file("stdin_bam", "-",
                                 pipe_from=f"{GOLDEN_DIR}/seqset/hiv_test.bam",
                                 data_format="bam")
        self.assertEqual(stats.paired_reads, 998)
        self.assertEqual(stats.unpaired_reads, 1)

    def test_indexed_bam(self):
        """Test indexed BAM file"""

        self.check(f"cp {GOLDEN_DIR}/spec/gatk/example_reads.bam {self.data_dir}/indexed_bam.bam")
        self.check(f"cp {GOLDEN_DIR}/spec/gatk/example_reads.bai {self.data_dir}/indexed_bam.bai")
        base_filename = f"{self.data_dir}/indexed_bam"
        bam_filename = base_filename + ".bam"
        idx_filename = base_filename + ".bai"
        self.assertEqual(self.index_filenames(base_filename), [idx_filename])
        stats = self.import_file("indexed_bam", bam_filename)
        self.assertEqual(stats.paired_reads, 32)
        self.assertEqual(stats.unpaired_reads, 1)

    def test_unindexed_cram(self):
        """Test unindexed CRAM file"""

        base_filename = f"{self.data_dir}/unindexed_cram"
        self.check(f"cp datasets/bams/e_coli/e_coli_test.cram {base_filename}.cram")
        self.assertEqual(self.index_filenames(base_filename), [],
                         msg="We're trying to test without indexes here")
        cram_filename = base_filename + ".cram"
        stats = self.import_file("unindexed_cram", cram_filename)
        self.assertEqual(stats.paired_reads, 53550)
        self.assertEqual(stats.unpaired_reads, 0)

    def test_indexed_cram(self):
        """Test indexed CRAM file"""

        base_filename = f"{self.data_dir}/indexed_cram"
        cram_filename = base_filename + ".cram"
        index_filename = cram_filename + ".crai"
        self.check(f"cp datasets/bams/e_coli/e_coli_test.cram {cram_filename}")

        # needed for samtools
        os.environ['PATH'] = os.environ['PATH'] + ":/share/software/bin"

        self.check(f"samtools index {cram_filename}")
        self.assertEqual(self.index_filenames(base_filename), [index_filename])
        stats = self.import_file("indexed_cram", cram_filename)
        self.assertEqual(stats.paired_reads, 53550)
        self.assertEqual(stats.unpaired_reads, 0)

    def test_multi_bam(self):
        """Test multiple bam files"""
        stats = self.import_file("multi_bam",
                                 f"{GOLDEN_DIR}/seqset/hiv_test.bam",
                                 f"{GOLDEN_DIR}/spec/gatk/example_reads.bam")
        self.assertEqual(stats.paired_reads,
                         998 # hiv_test
                         + 32) # example_reads
        self.assertEqual(stats.unpaired_reads,
                         1 # hiv_test
                         + 1) # example_reads

    def test_bam_and_fastq(self):
        """Test reading both a bam and fastq"""
        stats = self.import_file("multi_bam_fq",
                                 f"{GOLDEN_DIR}/seqset/hiv_test.bam",
                                 f"datasets/bams/pairing_test/pair_test1.fq")
        self.assertEqual(stats.paired_reads,
                         998 # hiv_test
                         + 0) # pair_test1.fq
        self.assertEqual(stats.unpaired_reads,
                         1 # hiv_test
                         + 25050) # pair_test1.fq

if __name__ == '__main__':
    unittest.main(verbosity=2)
