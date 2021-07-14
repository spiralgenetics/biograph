"""
biograph_test.py: Functional tests for biograph.
"""

from __future__ import print_function

import unittest
import os
import subprocess
from subprocess import (
    call,
    Popen,
    PIPE
)

import vcf


from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)

from python.functest.utils.fileops import (
    sha1_file
)

from python.functest.utils.defaults import (
    GOLDEN_DIR
)

from python.functest.utils.vcf import (
    complement
)

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class BioGraphTestCases(unittest.TestCase):
    """    unittest Test definitions follow """

    # keep pylint happy
    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit code matches (zero by default) """
        self.assertEqual(code, call(cmd, shell=True), "{0} did not exit {1}".format(cmd, code))

    def check_results(self, file1, file2):
        """ sha1 two files, throw if they don't match """
        self.assertEqual(sha1_file(test=self, file_name=file1), sha1_file(file2), "%s and %s do not match" % (file1, file2))

    # @unittest.skip(True)
    def test_00_basic(self):
        """
            Test essential functionality, including biograph --continue.
        """
        # self.cleanup = False
        # needed for samtools
        os.environ['PATH'] = os.environ['PATH'] + ":/share/software/bin"

        golden_fasta = '{0}/seqset/hiv_test.fa'.format(GOLDEN_DIR)
        golden_bam = '{0}/seqset/hiv_test.bam'.format(GOLDEN_DIR)
        generated_ref_dir = '{0}/hiv_test_ref'.format(self.data_dir)
        generated_biograph = '{0}/hiv_test.bg'.format(self.data_dir)
        generated_prefixes = '{0}/hiv_test.prefixes'.format(self.data_dir)
        generated_reads = '{0}/hiv_test.reads'.format(self.data_dir)

        # Basic test
        self.check("samtools view {0} |cut -f 10 |tr -d '\n' > {1}".format(golden_bam, generated_reads))
        self.check("bgbinary reference --in {0} --ref {1}".format(golden_fasta, generated_ref_dir))

        print('Starting biograph...')

        # Start a job
        pid = Popen([
            "bgbinary",
            "create",
            "--in", golden_bam,
            "--ref", generated_ref_dir,
            "--out", generated_biograph,
            "--tmp", self.data_dir + "/tmpdir",
            "--trim-after-portion", "1",
            "--max-corrections", "0",
            "--overrep-threshold", "0",
            "-f"
        ], shell=False, stdout=PIPE, stderr=PIPE)
        biograph_lines = pid.stderr.readlines()
        pid.stdout.close()
        pid.stderr.close()

        for expected in [
                '/seqset',
                '/qc/create_log.txt',
                '/qc/create_stats.json',
                '/qc/kmer_quality_report.html',
                '/metadata/bg_info.json'
        ]:
            thefile = generated_biograph + expected
            self.assertTrue(os.path.isfile(thefile), "{0} does not exist".format(thefile))

        for expected in [b"Importing reads\n",
                         b"Total reads imported: 999\n",
                         b"Running kmerization\n",
                         b"Correcting reads\n",
                         b"Generating BioGraph\n"]:
            self.assertIn(expected, biograph_lines, b"Can't find " + expected + b" in biograph output: " + b"".join(biograph_lines))
        self.assertRegex(biograph_lines[-2].decode(), "hiv_test.bg created\\.\n")

        self.check("bgbinary query --quiet --verbose --in {0} --query '' > {1}".format(generated_biograph + "/seqset", generated_prefixes))

        chrom = '>gi|9629357|ref|NC_001802.1| Human immunodeficiency virus 1, complete genome'
        reference_string = ""
        with open(golden_fasta, "r") as reference_stream:
            fasta_strings = reference_stream.readlines()
            self.assertEqual(fasta_strings[0], chrom + "\n", 'Has the tiny reference changed?  Expected "{0}" but found {1}'.format(chrom, fasta_strings[0]))
            reference_string = fasta_strings[1]
            self.assertEqual(len(reference_string), 9181, "Expected a 9181 base pair reference, but found {0}".format(len(reference_string)))

        reads = []
        with open(generated_reads, "r") as f:
            reads = f.readlines()

        with open(generated_prefixes, "r") as seqset_strings_stream:
            for line in seqset_strings_stream:
                line = line.strip()
                if line not in reads[0]:
                    self.assertIn(complement(line, True), reads[0]
                                  , 'Neither "{0}" nor its complement "{1}" from the SEQSET were found in the reads'.format(line, complement(line, True)))

    # @unittest.skip(True)
    def test_01_seqset_no_reads(self):
        """
            Test that an error is generated when all the reads are filtered out at anchoring
        """
        golden_bam = '{0}/seqset/tiny_test.bam'.format(GOLDEN_DIR)
        refdir = '/reference/e_coli_pairing_test/'
        generated_biograph = '{0}/tiny_test.bg'.format(self.data_dir)

        cmd = "bgbinary create --in {0} --ref {1} --out {2} --min-reads 0.6 --trim-after-portion 1 --max-corrections 0 --overrep-threshold 0 -f".format(golden_bam, refdir, generated_biograph)
        self.assertEqual(1, call(cmd, shell=True), "%s exited with zero but should not have done so" % cmd)

    # @unittest.skip(True)
    def test_02_merge(self):
        """
            Merge two biographs and test the merged result using the python bindings
        """
        # self.cleanup = False
        golden_fasta = '{0}/seqset/hiv_test.fa'.format(GOLDEN_DIR)
        golden_bam = '{0}/seqset/hiv_test.bam'.format(GOLDEN_DIR)
        generated_ref_dir = '{0}/hiv_test_ref'.format(self.data_dir)
        generated_biograph = '{0}/hiv_test.bg'.format(self.data_dir)
        generated_biograph2 = '{0}/hiv_test2.bg'.format(self.data_dir)
        merged_biograph = '{0}/hiv_merged.bg'.format(self.data_dir)
        create_stats_json = '{out}/qc/create_stats.json'.format(out=generated_biograph)
        create_stats_json2 = '{out}/qc/other_stats.json'.format(out=generated_biograph)

        # Generate two biographs
        self.check("bgbinary reference --in {0} --ref {1}".format(golden_fasta, generated_ref_dir))

        # Include a trailing / on --out to test for SEADEV-532
        self.check("bgbinary create --in {0} --ref {1} --out {2}/ --tmp {3} --trim-after-portion 1 --max-corrections 0 --overrep-threshold 0 --force"
                   .format(golden_bam, generated_ref_dir, generated_biograph, self.data_dir + "/tmpdir"))
        self.check("bgbinary create --in {0} --ref {1} --out {2} --tmp {3} --stats {4} --trim-after-portion 1 --max-corrections 0 --force --overrep-threshold 0 --id my_other_hiv_test"
                   .format(golden_bam, generated_ref_dir, generated_biograph2, self.data_dir + "/tmpdir", create_stats_json2))

        # stats JSON should exist
        self.assertTrue(os.path.exists(create_stats_json))
        self.assertTrue(os.path.exists(create_stats_json2))

        # Merge them
        self.check("bgbinary merge --out {0} --in {1} {2}".format(merged_biograph, generated_biograph, generated_biograph2))

        for expected in [
                '/seqset',
                '/qc/merge_log.txt',
                '/qc/merge_stats.json',
                '/metadata/bg_info.json'
        ]:
            thefile = merged_biograph + expected
            self.assertTrue(os.path.isfile(thefile), "{0} does not exist".format(thefile))

        # Verify the reference
        from biograph import BioGraph, Reference

        ref = Reference(generated_ref_dir)

        self.assertEqual(len(ref.scaffolds), 1)
        self.assertEqual(len(ref.supercontigs), 1)

        self.assertEqual(ref.supercontigs[0].size, 9180)
        self.assertEqual(ref.scaffolds[0], 'gi|9629357|ref|NC_001802.1|')

        with self.assertRaises(RuntimeError):
            ref.make_range(ref.scaffolds[0], 1, 9183)

        # Should not raise
        ref.make_range(ref.scaffolds[0], 0, 9181)

        # Verify the merged biograph
        my_bg = BioGraph(merged_biograph)

        self.assertEqual(my_bg.metadata.accession_id, 'hiv_test+my_other_hiv_test')

        # both samples should be identical
        self.assertEqual(len(my_bg.metadata.samples), 2)

        # metadata update
        self.check("bgbinary metadata --in {0} --accession-id different --sample-id hiv_test=foo".format(merged_biograph))
        my_new_bg = BioGraph(merged_biograph)
        self.assertEqual(my_new_bg.metadata.accession_id, 'different')
        self.assertEqual(my_bg.metadata.samples['hiv_test'], my_new_bg.metadata.samples['foo'])

    # @unittest.skip(True)
    def test_03_seqset_python_bindings(self):
        """
            Check the python bindings on a bigger reference
        """
        ref = '/reference/human_g1k_v37'

        #self.cleanup = False

        from biograph import Reference
        ref = Reference(ref)

        # 84 scaffolds, 338 supercontigs in human_g1k_v37
        self.assertEqual(len(ref.scaffolds), 84)
        self.assertEqual(len(ref.supercontigs), 338)

        refrange = ref.make_range("2", 100000, 200003, False)
        self.assertEqual(refrange.start, 100000)
        self.assertEqual(refrange.end, 200003)
        self.assertEqual(refrange.size, 100003)
        self.assertEqual(refrange.scaffold, "2")

        with self.assertRaises(RuntimeError):
            refrange = ref.make_range("20", 1, 177954, True)

        # should NOT raise
        refrange = ref.make_range("20", 1, 177954, False)

        refrange = ref.make_range("20", 60127, 60137, True)
        self.assertEqual(str(refrange.sequence), 'ACCATGGACC')

    # @unittest.skip(True)
    def test_paired_fastq(self):
        """
            Test FASTQ read pairing
        """
        # self.cleanup = False

        refdir = '/reference/e_coli_pairing_test/'
        fastq1 = '/share/bams/pairing_test/pair_test1.fq'
        fastq2 = '/share/bams/pairing_test/pair_test2.fq'
        test_sequences = '{out}/test.sequences'.format(out=self.data_dir)
        test_biograph = '{out}/pair_test.bg'.format(out=self.data_dir)
        golden = '/share/bams/pairing_test/20170804/pair_test.bg/seqset'
        golden_sequences = '{out}/golden.sequences'.format(out=self.data_dir)

        self.check("bgbinary create --in {0} --pair {1} --ref {2} --out {3} --trim-after-portion 1 --max-corrections 0 --overrep-threshold 0".format(fastq1, fastq2, refdir, test_biograph))
        self.check("bgbinary query --in {0} --query '' --verbose > {1} 2>/dev/null".format(golden, golden_sequences))
        self.check("bgbinary query --in {0} --query '' --verbose > {1} 2>/dev/null".format(test_biograph + "/seqset", test_sequences))

        self.assertEqual(sha1_file(test=self, file_name=golden_sequences), sha1_file(test_sequences), "%s and %s do not match" % (golden_sequences, test_sequences))

        # Can we find the SNP?
        from biograph import BioGraph

        my_bg = BioGraph(test_biograph)

        self.assertEqual(my_bg.seqset.size(), 196772)
        self.assertEqual("CCACACACGACGGTTACGGTTTTCTCCCATTGCCAGCAGACCACGAATGGCAATGGACAGCCACATGGAGATCACCGGCTTCAGGGAGGCA", str(my_bg.seqset.get_entry_by_id(65590).sequence()))

    # @unittest.skip(True)
    def test_paired_interleaved_fastq(self):
        """
            Test FASTQ interleaved read pairing
        """
        # self.cleanup = False

        refdir = '/reference/e_coli_pairing_test/'
        fastq = '/share/bams/pairing_test/interleaved.fq'
        test_sequences = '{out}/test.sequences'.format(out=self.data_dir)
        test_biograph = '{out}/pair_test.bg'.format(out=self.data_dir)
        golden = '/share/bams/pairing_test/20170804/pair_test.bg/seqset'
        golden_sequences = '{out}/golden.sequences'.format(out=self.data_dir)

        self.check("bgbinary create --in {0} --ref {1} --out {2} --interleaved --trim-after-portion 1 --max-corrections 0 --overrep-threshold 0".format(fastq, refdir, test_biograph))
        self.check("bgbinary query --in {0} --query '' --verbose > {1} 2>/dev/null".format(golden, golden_sequences))
        self.check("bgbinary query --in {0} --query '' --verbose > {1} 2>/dev/null".format(test_biograph + "/seqset", test_sequences))

        self.assertEqual(sha1_file(test=self, file_name=golden_sequences), sha1_file(test_sequences), "%s and %s do not match" % (golden_sequences, test_sequences))

        # Can we find the SNP?
        from biograph import BioGraph

        my_bg = BioGraph(test_biograph)

        self.assertEqual(my_bg.seqset.size(), 196772)
        self.assertEqual("CCACACACGACGGTTACGGTTTTCTCCCATTGCCAGCAGACCACGAATGGCAATGGACAGCCACATGGAGATCACCGGCTTCAGGGAGGCA", str(my_bg.seqset.get_entry_by_id(65590).sequence()))

    # @unittest.skip(True)
    def test_variants_vcf(self):
        """
            Test variants VCF generation
        """
        # self.cleanup = False

        refdir = '/reference/e_coli_pairing_test/'
        fastq1 = '/share/bams/pairing_test/pair_test1.fq'
        fastq2 = '/share/bams/pairing_test/pair_test2.fq'
        test_biograph = '{out}/pair_test.bg'.format(out=self.data_dir)
        test_vcf = '{out}/variants.vcf'.format(out=self.data_dir)
        test_vcf_with_AID = '{out}/variants-assemblies.vcf'.format(out=self.data_dir)
        test_assemblies = '{out}/assemblies.csv'.format(out=self.data_dir)
        create_stats_json = '{out}/qc/create_stats.json'.format(out=test_biograph)
        variants_stats_json = '{out}/qc/variants_stats.json'.format(out=test_biograph)
        variants_stats_json2 = '{out}/qc/some_other_variants_stats.json'.format(out=test_biograph)

        self.check("bgbinary create --in {0} --pair {1} --ref {2} --out {3}".format(fastq1, fastq2, refdir, test_biograph))

        # stats JSON should exist
        self.assertTrue(os.path.exists(create_stats_json))

        # With and without assemblies
        self.check("bgbinary discovery --rvg-exclude=false --in {0} --ref {1} --out {2} --enable-pop-tracer=false".format(test_biograph, refdir, test_vcf))
        self.check("bgbinary discovery --rvg-exclude=false --in {0} --ref {1} --out {2} --stats {3} --assemblies-out {4} --enable-pop-tracer=false -f".format(test_biograph, refdir, test_vcf_with_AID, variants_stats_json2, test_assemblies))

        # stats JSONs should exist
        self.assertTrue(os.path.exists(variants_stats_json))
        self.assertTrue(os.path.exists(variants_stats_json2))

        # Count header lines without output assemblies
        self.assertEqual(subprocess.check_output("grep -c ^# {0}".format(test_vcf), shell=True).strip(), b"29")
        # no AID header line expected
        self.check("grep '^##INFO=<ID=AID' {0}".format(test_vcf), 1)
        # No variants with AID in the INFO
        self.check("grep -c AID {0}".format(test_vcf), 1)

        # Count header lines with output assemblies
        self.assertEqual(subprocess.check_output("grep -c ^# {0}".format(test_vcf_with_AID), shell=True).strip(), b"30")
        # explicit AID header
        self.check("grep -q '^##INFO=<ID=AID' {0}".format(test_vcf_with_AID))
        # 6 variants + header with AID
        self.assertEqual(subprocess.check_output("grep -c AID {0}".format(test_vcf_with_AID), shell=True).strip(), b"19")

        for the_vcf in (test_vcf, test_vcf_with_AID):
            # reference should match
            self.check("grep '^##reference={0}$' {1}".format(refdir, the_vcf))

            # Every one of these header lines must exist
            for header in ('##fileformat', '##fileDate', '##source', '##INFO', '##FORMAT', '##ALT', '##contig', '#CHROM'):
                self.check("grep '^{0}' {1}".format(header, the_vcf))

            # Should be valid VCF
            with open(the_vcf) as vcf_file:
                my_vcf = vcf.Reader(vcf_file)
                cnt = 0
                for _ in my_vcf:
                    cnt += 1
                self.assertEqual(cnt, 18)

    # @unittest.skip(True)
    def test_variants_vcf_multithreaded(self):
        """
            Ensure multithreaded runs always have the same number of variants when the pop tracer is enabled.
        """
        refdir = '/reference/e_coli_k12_ASM584v1/'
        bam = '/share/bams/e_coli/e_coli_test.bam'
        test_biograph = '{out}/e_coli_test.bg'.format(out=self.data_dir)
        test_vcf = '{out}/variants.vcf'.format(out=self.data_dir)

        self.check("bgbinary create --in {0} --ref {1} --out {2}".format(bam, refdir, test_biograph))

        for i in range(10):
            self.check("bgbinary discovery --rvg-exclude=false --in {0} --ref {1} --out {2} -f".format(test_biograph, refdir, test_vcf))

            # Count header lines
            self.assertEqual(subprocess.check_output("grep -c ^# {0}".format(test_vcf), shell=True).strip(), b"31", "Incorrect number of header lines at iteration #{0}".format(i + 1))

            # Count variants
            self.assertEqual(subprocess.check_output("grep -cv ^# {0}".format(test_vcf), shell=True).strip(), b"140", "Incorrect number of vcf entries at iteration #{0}".format(i + 1))

if __name__ == '__main__':
    unittest.main(verbosity=2)
