"""
spec.py: Functional tests for spec.
"""

from __future__ import print_function

import unittest
import subprocess
import os
import sys
import random
from time import time

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

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class SpecTestCases(unittest.TestCase):
    """    unittest Test definitions follow """

    # keep pylint happy
    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def check(self, cmd):
        """ Run a shell command, assert exit is zero. """
        self.assertEqual(0, subprocess.call(cmd, shell=True), "%s exited non-zero" % cmd)

    def check_results(self, file1, file2):
        """ assert two files sha the same """
        self.assertEqual(sha1_file(test=self, file_name=file1), sha1_file(file2), "%s and %s do not match" % (file1, file2))

    def test_szip(self):
        """
            szip test
        """

        # self.cleanup = False
        golden_file = '%s/quick_multisize.fq' % GOLDEN_DIR
        compressed_file = '%s/qm.sz' % self.data_dir
        uncompressed_file = '%s/qm.fq' % self.data_dir

        # Basic test
        self.check("szip %s %s" % (golden_file, compressed_file))
        self.check("szip -d %s %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)
        self.assertTrue(os.path.getsize(compressed_file) < os.path.getsize(uncompressed_file), "Compressed file is larger than uncompressed file.")

        # Refuse to overwrite output without -f
        self.assertEqual(1, subprocess.call("szip -d %s %s 2>/dev/null" % (compressed_file, uncompressed_file), shell=True), "szip should refuse to overwrite output without -f but didn't")
        self.assertEqual(0, subprocess.call("szip -df %s %s" % (compressed_file, uncompressed_file), shell=True), "szip should overwrite output with -f but didn't")
        self.check_results(uncompressed_file, golden_file)

        # Compressed results should match
        self.check("szip %s %s_2" % (golden_file, compressed_file))
        self.check_results(compressed_file, "%s_2" % compressed_file)

        # STDOUT
        self.check("szip %s > %s" % (golden_file, compressed_file))
        self.check("szip -d %s > %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

        # STDIN and STDOUT
        self.check("szip < %s > %s" % (golden_file, compressed_file))
        self.check("szip -d < %s > %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

        # STDIN - and STDOUT -
        self.check("szip - - < %s > %s" % (golden_file, compressed_file))
        self.check("szip -d - - < %s > %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

        # STDIN - and STDOUT -
        self.check("szip --in - --out - < %s > %s" % (golden_file, compressed_file))
        self.check("szip -d - --out - < %s > %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

        # STDIN | and STDOUT
        self.check("cat %s | szip > %s" % (golden_file, compressed_file))
        self.check("cat %s | szip -d > %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

        # Let's mix it up
        self.check("cat %s | szip | tee %s > /dev/null" % (golden_file, compressed_file))
        self.check("szip -df %s %s" % (compressed_file, uncompressed_file))
        self.check_results(uncompressed_file, golden_file)

    #pylint: disable=too-many-locals, too-many-statements
    def test_00_spec(self):
        """
            spec2bam and bam2spec
        """

        # self.cleanup = False
        golden_bam = '%s/spec/test.bam' % GOLDEN_DIR
        golden_spec = '%s/spec/test.spec' % GOLDEN_DIR
        golden_stats = '%s/spec/stats.txt' % GOLDEN_DIR
        ref = '/share/flat_refs/hg19.flat'
        wrong_ref = '/share/flat_refs/e_coli_k12_ASM584v1.ref'
        test_bam = '%s/test.bam' % self.data_dir
        test_spec = '%s/test.spec' % self.data_dir
        uncompressed_bam = '%s/uncompressed.bam' % self.data_dir
        binned_spec = '%s/test-binned.spec' % self.data_dir
        binned_bam = '%s/test-binned.bam' % self.data_dir
        discard_spec = '%s/test-discard.spec' % self.data_dir
        discard_bam = '%s/test-discard.bam' % self.data_dir
        keep_spec = '%s/test-keep.spec' % self.data_dir
        keep_bam = '%s/test-keep.bam' % self.data_dir
        stats = '%s/stats.txt' % self.data_dir

        # Basic test
        self.check("bam2spec --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))
        self.check("spec2bam --in %s --ref %s --out %s" % (test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)

        # Invalid inputs
        self.assertEqual(1, subprocess.call("bam2spec --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "bam2spec_integrate should refuse to encode spec files")
        self.assertEqual(1, subprocess.call("spec2bam --in %s --ref %s --out %s 2>/dev/null" % (test_bam, ref, test_spec), shell=True), "spec2bam should refuse to encode bam files")

        self.assertEqual(1, subprocess.call("bam2spec --threads foo --in %s --ref %s --out %s 2>/dev/null" % (test_bam, ref, test_spec), shell=True), "bam2spec should invalid thread count")
        self.assertEqual(1, subprocess.call("spec2bam --threads 65535 --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "spec2bam should refuse huge thread count")

        # Bad reference
        self.assertEqual(1, subprocess.call("bam2spec --in %s --ref %s --out %s 2>/dev/null" % (test_bam, wrong_ref, test_spec), shell=True), "bam2spec should refuse to use the wrong reference")

        # Force overwrite
        self.assertEqual(1, subprocess.call("bam2spec --in %s --ref %s --out %s 2>/dev/null" % (golden_bam, ref, test_spec), shell=True), "bam2spec should refuse to overwrite existing spec file")
        self.check("bam2spec --force --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))

        self.assertEqual(1, subprocess.call("spec2bam --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "spec2bam should refuse to overwrite existing spec file")
        self.check("spec2bam --force --in %s --ref %s --out %s" % (test_spec, ref, test_bam))

        # Above tests should not change previous test results
        self.check_results(golden_bam, test_bam)

        # Uncompressed bam
        self.check("spec2bam --uncompressed --in %s --ref %s --out %s" % (test_spec, ref, uncompressed_bam))
        self.assertTrue(os.path.getsize(uncompressed_bam) > os.path.getsize(golden_bam), "Compressed file is larger than uncompressed file.")

        # Quality binning
        self.check("bam2spec --bin-qualities --in %s --ref %s --out %s" % (golden_bam, ref, binned_spec))
        self.assertTrue(os.path.getsize(binned_spec) < os.path.getsize(golden_spec), "Binned spec file is larger than uncompressed file.")

        self.check("spec2bam --in %s --ref %s --out %s" % (binned_spec, ref, binned_bam))
        self.assertTrue(os.path.getsize(binned_bam) < os.path.getsize(golden_bam), "Binned bam file is larger than uncompressed file.")

        # Discard some attributes
        self.check("bam2spec --in %s --ref %s --out %s --discard-tags MD" % (golden_bam, ref, discard_spec))
        self.assertTrue(os.path.getsize(discard_spec) < os.path.getsize(golden_spec), "Discard spec file is larger than golden spec.")

        self.check("spec2bam --in %s --ref %s --out %s" % (discard_spec, ref, discard_bam))
        self.assertTrue(os.path.getsize(discard_bam) < os.path.getsize(golden_bam), "Discard bam file is larger than golden bam.")

        # Discard most attributes
        self.check("bam2spec --in %s --ref %s --out %s --keep-tags X0:AM" % (golden_bam, ref, keep_spec))
        self.assertTrue(os.path.getsize(keep_spec) < os.path.getsize(golden_spec), "Keep spec file is larger than golden spec.")
        self.assertTrue(os.path.getsize(keep_spec) < os.path.getsize(discard_spec), "Keep spec file is larger than discard spec.")

        self.check("spec2bam --in %s --ref %s --out %s" % (keep_spec, ref, keep_bam))
        self.assertTrue(os.path.getsize(keep_bam) < os.path.getsize(golden_bam), "Keep bam file is larger than golden bam.")
        self.assertTrue(os.path.getsize(keep_bam) < os.path.getsize(discard_bam), "Keep bam file is larger than discard bam.")

        # Discard all attributes
        self.check("bam2spec --force --discard-tags all --in %s --ref %s --out %s" % (golden_bam, ref, discard_spec))
        self.assertTrue(os.path.getsize(discard_spec) < os.path.getsize(golden_spec), "Discard all spec file is larger than golden spec.")

        self.check("bam2spec -f --keep-tags none --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))
        self.assertTrue(os.path.getsize(test_spec) < os.path.getsize(golden_spec), "Keep none spec file is larger than golden spec.")

        self.check("spec2bam -f --in %s --ref %s --out %s" % (discard_spec, ref, discard_bam))
        self.assertTrue(os.path.getsize(discard_bam) < os.path.getsize(golden_bam), "Discard all bam file is larger than golden bam.")

        self.check("spec2bam --force --in %s --ref %s --out %s" % (discard_spec, ref, test_bam))
        self.assertTrue(os.path.getsize(test_bam) < os.path.getsize(golden_bam), "Keep none bam file is larger than golden bam.")

        self.check_results(discard_bam, test_bam)

        # Stats
        self.check("bam2spec -f --stats --in %s --ref %s --out %s > %s" % (golden_bam, ref, test_spec, stats))
        self.check_results(golden_stats, stats)

        # Encryption
        self.check("SPEC_KEY=foo bam2spec --force --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))

        self.check("spec2bam --key foo --force --in %s --ref %s --out %s" % (test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)

        # Command line overrides environment variable
        self.check("SPEC_KEY=bar spec2bam --key foo --force --in %s --ref %s --out %s" % (test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)
        self.assertEqual(2, subprocess.call("spec2bam --key bar --force --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "spec2bam should fail to decrypt with the wrong key")

        # SPEC_REF
        self.check("SPEC_REF=%s bam2spec --force --in %s --out %s" % (ref, golden_bam, test_spec))
        self.check("SPEC_REF=%s spec2bam --force --in %s --out %s" % (ref, test_spec, test_bam))
        self.check_results(golden_bam, test_bam)

        # Command line overrides environment variable
        self.assertEqual(2, subprocess.call("SPEC_REF=%s spec2bam --force --in %s --ref foo --out %s 2>/dev/null" % (ref, test_spec, test_bam), shell=True), "spec2bam should fail to decode with wrong ref")

        # write-bai
        self.check("bam2spec --in %s --ref %s --out %s --write-bai --force" % (golden_bam, ref, test_spec))
        self.assertTrue(os.path.exists("%s.bai" % test_spec), "--write-bai did not create a .spec.bai file")
        self.assertTrue(os.path.getsize("%s.bai" % test_spec) > 0, "%s.bai is empty" % test_spec)

    def test_spec_decoder(self):
        """
            spec_decoder: should only run as spec2bam
        """

        # self.cleanup = False
        golden_bam = '%s/spec/test.bam' % GOLDEN_DIR
        golden_spec = '%s/spec/test.spec' % GOLDEN_DIR
        golden_stats = '%s/spec/stats.txt' % GOLDEN_DIR
        ref = '/share/flat_refs/hg19.flat'
        golden_fasta = '%s/spec/gatk/exampleFASTA.fasta' % GOLDEN_DIR
        test_bam = '%s/test.bam' % self.data_dir
        test_spec = '%s/test.spec' % self.data_dir
        uncompressed_bam = '%s/uncompressed.bam' % self.data_dir
        binned_spec = '%s/test-binned.spec' % self.data_dir
        binned_bam = '%s/test-binned.bam' % self.data_dir
        discard_spec = '%s/test-discard.spec' % self.data_dir
        discard_bam = '%s/test-discard.bam' % self.data_dir
        keep_spec = '%s/test-keep.spec' % self.data_dir
        keep_bam = '%s/test-keep.bam' % self.data_dir
        stats = '%s/stats.txt' % self.data_dir

        decoder = os.getcwd() + '/modules/spec/spec2bam'

        # Basic test
        self.check("bam2spec --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))
        self.check("%s --in %s --ref %s --out %s" % (decoder, test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)

        # Invalid inputs
        self.assertEqual(1, subprocess.call("bam2spec --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "bam2spec should refuse to encode spec files")
        self.assertEqual(1, subprocess.call("%s --in %s --ref %s --out %s 2>/dev/null" % (decoder, test_bam, ref, test_spec), shell=True), "spec2bam should refuse to encode bam files")

        self.assertEqual(1, subprocess.call("bam2spec --threads foo --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "bam2spec should invalid thread count")
        self.assertEqual(1, subprocess.call("%s --threads 65535 --in %s --ref %s --out %s 2>/dev/null" % (decoder, test_bam, ref, test_spec), shell=True), "spec2bam should refuse huge thread count")

        # Force overwrite
        self.assertEqual(1, subprocess.call("bam2spec --in %s --ref %s --out %s 2>/dev/null" % (golden_bam, ref, test_spec), shell=True), "bam2spec should refuse to overwrite existing spec file")
        self.check("bam2spec --force --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))

        self.assertEqual(1, subprocess.call("%s --in %s --ref %s --out %s 2>/dev/null" % (decoder, test_spec, ref, test_bam), shell=True), "spec2bam should refuse to overwrite existing spec file")
        self.check("%s --force --in %s --ref %s --out %s" % (decoder, test_spec, ref, test_bam))

        # Above tests should not change previous test results
        self.check_results(golden_bam, test_bam)

        # Uncompressed bam
        self.check("%s --uncompressed --in %s --ref %s --out %s" % (decoder, test_spec, ref, uncompressed_bam))
        self.assertTrue(os.path.getsize(uncompressed_bam) > os.path.getsize(golden_bam), "Compressed file is larger than uncompressed file.")

        # Quality binning
        self.check("bam2spec --bin-qualities --in %s --ref %s --out %s" % (golden_bam, ref, binned_spec))
        self.assertTrue(os.path.getsize(binned_spec) < os.path.getsize(golden_spec), "Binned spec file is larger than uncompressed file.")

        self.check("%s --in %s --ref %s --out %s" % (decoder, binned_spec, ref, binned_bam))
        self.assertTrue(os.path.getsize(binned_bam) < os.path.getsize(golden_bam), "Binned bam file is larger than uncompressed file.")

        # Discard some attributes
        self.check("bam2spec --in %s --ref %s --out %s --discard-tags MD" % (golden_bam, ref, discard_spec))
        self.assertTrue(os.path.getsize(discard_spec) < os.path.getsize(golden_spec), "Discard spec file is larger than golden spec.")

        self.check("%s --in %s --ref %s --out %s" % (decoder, discard_spec, ref, discard_bam))
        self.assertTrue(os.path.getsize(discard_bam) < os.path.getsize(golden_bam), "Discard bam file is larger than golden bam.")

        # Discard most attributes
        self.check("bam2spec --in %s --ref %s --out %s --keep-tags X0:AM" % (golden_bam, ref, keep_spec))
        self.assertTrue(os.path.getsize(keep_spec) < os.path.getsize(golden_spec), "Keep spec file is larger than golden spec.")
        self.assertTrue(os.path.getsize(keep_spec) < os.path.getsize(discard_spec), "Keep spec file is larger than discard spec.")

        self.check("%s --in %s --ref %s --out %s" % (decoder, keep_spec, ref, keep_bam))
        self.assertTrue(os.path.getsize(keep_bam) < os.path.getsize(golden_bam), "Keep bam file is larger than golden bam.")
        self.assertTrue(os.path.getsize(keep_bam) < os.path.getsize(discard_bam), "Keep bam file is larger than discard bam.")

        # Discard all attributes
        self.check("bam2spec --force --discard-tags all --in %s --ref %s --out %s" % (golden_bam, ref, discard_spec))
        self.assertTrue(os.path.getsize(discard_spec) < os.path.getsize(golden_spec), "Discard all spec file is larger than golden spec.")

        self.check("bam2spec -f --keep-tags none --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))
        self.assertTrue(os.path.getsize(test_spec) < os.path.getsize(golden_spec), "Keep none spec file is larger than golden spec.")

        self.check("%s -f --in %s --ref %s --out %s" % (decoder, discard_spec, ref, discard_bam))
        self.assertTrue(os.path.getsize(discard_bam) < os.path.getsize(golden_bam), "Discard all bam file is larger than golden bam.")

        self.check("%s --force --in %s --ref %s --out %s" % (decoder, discard_spec, ref, test_bam))
        self.assertTrue(os.path.getsize(test_bam) < os.path.getsize(golden_bam), "Keep none bam file is larger than golden bam.")

        self.check_results(discard_bam, test_bam)

        # Stats
        self.check("bam2spec -f --stats --in %s --ref %s --out %s > %s" % (golden_bam, ref, test_spec, stats))
        self.check_results(golden_stats, stats)

        # Encryption
        self.check("SPEC_KEY=foo bam2spec --force --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))

        self.check("spec2bam --key foo --force --in %s --ref %s --out %s" % (test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)

        # Command line overrides environment variable
        self.check("SPEC_KEY=bad spec2bam --key foo --force --in %s --ref %s --out %s" % (test_spec, ref, test_bam))
        self.check_results(golden_bam, test_bam)

        self.assertEqual(2, subprocess.call("spec2bam --key bar --force --in %s --ref %s --out %s 2>/dev/null" % (test_spec, ref, test_bam), shell=True), "spec2bam should fail to decrypt with the wrong key")

        # No other binaries should run
        self.check("ln -s %s %s/bam2spec" % (decoder, self.data_dir))
        self.check("ln -s %s %s/fasta2ref" % (decoder, self.data_dir))

        self.assertEqual(1, subprocess.call("%s/bam2spec --in %s --ref %s --out %s/does_not_exist.spec 2>/dev/null" % (self.data_dir, golden_bam, ref, self.data_dir), shell=True), "bam2spec must not work with decoder binary")
        self.assertFalse(os.path.isfile('%s/does_not_exist.spec' % self.data_dir), "decoder-only bam2spec created a spec file!")
        self.assertEqual(1, subprocess.call("%s/fasta2ref %s %s/out.ref 2>/dev/null" % (self.data_dir, golden_fasta, self.data_dir), shell=True), "fasta2ref must not work with decoder binary")

    @unittest.skipIf(not os.path.isfile('/opt/spiral/lib/spec-1.120.jar'), 'specstream not found')
    def test_zz_picard_read(self):
        """
            spec support for Picard
        """
        def ghetto_fastq_sort(infile, outfile):
            """ sort a fastq with paste and tr """
            self.check('cat %s | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > %s' % (infile, outfile))

        # self.cleanup = False

        ref = '/share/flat_refs/hg19.flat'

        golden_bam = '%s/spec/test.bam' % GOLDEN_DIR
        golden_sam = '%s/spec/test.sam' % GOLDEN_DIR
        golden_spec = '%s/spec/test.spec' % GOLDEN_DIR

        test_bam = '%s/test.bam' % self.data_dir
        test_sam = '%s/test.sam' % self.data_dir
        test_sam2 = '%s/test2.sam' % self.data_dir
        test_spec = '%s/test.spec' % self.data_dir
        test_fastq = '%s/test.fastq' % self.data_dir
        test_fastq2 = '%s/test2.fastq' % self.data_dir
        test_stats = '%s/test.txt' % self.data_dir
        test_stats2 = '%s/test2.txt' % self.data_dir

        java_home = "/usr/lib/jvm/java-1.7.0-openjdk-amd64"
        java = "JAVA_HOME=%s %s/bin/java" % (java_home, java_home)

        picard = 'SPEC_REF="%s" %s -jar /opt/spiral/lib/picard-1.130.jar' % (ref, java)

        # ViewSam
        self.check("%s ViewSam I=%s RECORDS_ONLY=true > %s 2>/dev/null" % (picard, golden_spec, test_sam))
        self.check("%s ViewSam I=%s RECORDS_ONLY=true > %s 2>/dev/null" % (picard, golden_bam, test_sam2))
        self.check_results(test_sam, test_sam2)
        self.check_results(test_sam, golden_sam)

        # SamFormatConverter
        self.check("%s SamFormatConverter I=%s O=%s 2>/dev/null" % (picard, golden_spec, test_bam))
        self.check("%s ViewSam I=%s RECORDS_ONLY=true > %s 2>/dev/null" % (picard, test_bam, test_sam))
        self.check_results(test_sam, test_sam2)

        # CollectQualityYieldMetrics
        self.check("%s CollectQualityYieldMetrics I=%s O=%s 2>/dev/null" % (picard, golden_spec, test_stats))
        self.check("%s CollectQualityYieldMetrics I=%s O=%s 2>/dev/null" % (picard, golden_bam, test_stats2))
        self.check("diff -I '^#' %s %s" % (test_stats, test_stats2))

        # SamToFastq
        self.check("%s SamToFastq I=%s F=%s.tmp F2=%s2.tmp 2>/dev/null" % (picard, golden_spec, test_fastq, test_fastq))
        ghetto_fastq_sort(test_fastq + ".tmp", test_fastq)
        self.check("%s SamToFastq I=%s F=%s.tmp F2=%s2.tmp 2>/dev/null" % (picard, golden_bam, test_fastq2, test_fastq2))
        ghetto_fastq_sort(test_fastq2 + ".tmp", test_fastq2)
        self.check_results(test_fastq, test_fastq2)

        # Encryption
        test_spec = "%s/encrypted.spec" % self.data_dir
        self.check("SPEC_KEY=baz bam2spec --in %s --out %s --ref %s  --no-match-reference > /dev/null 2>&1" % (golden_bam, test_spec, ref))
        self.assertEqual(1, subprocess.call("%s ViewSam I=%s RECORDS_ONLY=true > %s 2>/dev/null" % (picard, test_spec, test_sam2), shell=True), "%s appears to be unencrypted" % test_spec)

        self.check("SPEC_KEY=baz %s ViewSam I=%s RECORDS_ONLY=true > %s 2>/dev/null" % (picard, golden_bam, test_sam2))
        self.check_results(test_sam, test_sam2)

    @unittest.skipIf(not os.path.isfile('/opt/spiral/lib/spec-1.120.jar'), 'specstream not found')
    def test_zz_picard_write(self):
        """
            spec write support for Picard
        """
        # self.cleanup = False
        fastq = '/src/golden/ftest/E_coli_phred64.fq'
        ref = '/share/flat_refs/e_coli_k12_ASM584v1.ref'
        golden_sam = '%s/spec/e_coli.sam' % GOLDEN_DIR

        java_home = "/usr/lib/jvm/java-1.7.0-openjdk-amd64"
        java = "JAVA_HOME=%s %s/bin/java" % (java_home, java_home)
        picard = 'SPEC_REF="%s" %s -jar /opt/spiral/lib/picard-1.130.jar' % (ref, java)

        test_pat = '%s/test' % self.data_dir

        self.check("%s FastqToSam FASTQ=%s O=%s1.spec SAMPLE_NAME='e.coli' 2>/dev/null" % (picard, fastq, test_pat))
        self.check("spec2bam --in %s1.spec --out %s1.bam --ref %s" % (test_pat, test_pat, ref))
        self.check("samtools view -h %s1.bam > %s1.sam" % (test_pat, test_pat))
        self.check_results("%s1.sam" % test_pat, golden_sam)

        # encryption
        self.check("SPEC_KEY=bar %s FastqToSam FASTQ=%s O=%s2.spec SAMPLE_NAME='e.coli' 2>/dev/null" % (picard, fastq, test_pat))

        # fail, wrong key
        self.assertEqual(2, subprocess.call("spec2bam --in %s2.spec --out %s2.bam --ref %s > /dev/null 2>&1" % (test_pat, test_pat, ref), shell=True), "%s.spec does not appear to be encrypted." % test_pat)

        self.check("SPEC_KEY=bar spec2bam --in %s2.spec --out %s2.bam --ref %s" % (test_pat, test_pat, ref))
        self.check("samtools view -h %s2.bam > %s2.sam" % (test_pat, test_pat))
        self.check_results("%s2.sam" % test_pat, golden_sam)

    def test_index(self):
        """
            Validate that indexing works
        """

        # Set all the variables
        # self.cleanup = False
        golden_bam = '/share/bams/NA12878_S1-sample.bam'
        ref = '/share/flat_refs/hg19.flat'
        test_header = '%s/header.txt' % self.data_dir
        test_spec = '%s/test.spec' % self.data_dir
        test_bam = '%s/test.bam' % self.data_dir
        test_sam1 = '%s/test1.sam' % self.data_dir
        test_sam2 = '%s/test2.sam' % self.data_dir

        # Extract header
        print("Extracting header")
        self.check("samtools view -H %s > %s" % (golden_bam, test_header))
        chrs = []
        chr_lens = {}
        with open(test_header, "r") as ins:
            for line in ins:
                by_tabs = line.split('\t')
                if by_tabs[0] != '@SQ':
                    continue
                by_colon = by_tabs[1].split(':')
                the_chr = by_colon[1]
                by_colon = by_tabs[2].split(':')
                the_len = int(by_colon[1])
                chrs += [the_chr]
                chr_lens[the_chr] = the_len

        # Make indexed spec, removing any old version
        print("Making spec file")
        self.check("bam2spec --in %s --ref %s --out %s" % (golden_bam, ref, test_spec))

        for i in range(10):
            # Pick a random chromosome + range
            region_chr = random.choice(chrs)
            tot_len = chr_lens[region_chr]
            the_len = random.randrange(1, int((tot_len - 300) / random.randrange(3, 100)) + 2)
            region_start = random.randrange(1, tot_len - the_len - 10)
            region_end = region_start + the_len

            region_str = "%s:%d-%d" % (region_chr, region_start, region_end)
            print("%d: Checking %s" % (i, region_str))
            self.check("spec2bam -f --in %s --ref %s --out %s --region %s" % (test_spec, ref, test_bam, region_str))
            self.check("samtools view %s > %s" % (test_bam, test_sam1))
            self.check("samtools view %s %s > %s" % (golden_bam, region_str, test_sam2))
            self.check_results(test_sam1, test_sam2)

    @unittest.skipIf(not os.path.isfile('/opt/spiral/lib/spec-1.120.jar'), 'specstream not found')
    def test_zz_gatk_read(self):
        """
            spec read support for GATK
        """
        # self.cleanup = False

        ref = '/share/flat_refs/hg19.flat'

        golden_spec = '%s/spec/gatk/example.spec' % GOLDEN_DIR
        golden_sam = '%s/spec/gatk/example_reads.sam' % GOLDEN_DIR
        golden_bam = '%s/spec/gatk/example_reads.bam' % GOLDEN_DIR
        golden_loci_bam = '%s/spec/gatk/loci.bam' % GOLDEN_DIR
        golden_clipped_bam = '%s/spec/gatk/example_reads_clipped.bam' % GOLDEN_DIR
        golden_fasta = '%s/spec/gatk/exampleFASTA.fasta' % GOLDEN_DIR

        java_home = "/usr/lib/jvm/java-1.7.0-openjdk-amd64"
        java = "JAVA_HOME=%s %s/bin/java" % (java_home, java_home)
        gatk = "/share/gatk_resources/gatk-spiral-3.3-0.jar"

        test = 0
        for jar in ["/opt/spiral/lib/spec-1.120.jar", "/opt/spiral/lib/spec-1.130.jar"]:

            cmdline = "%s -cp %s:%s org.broadinstitute.gatk.engine.CommandLineGATK --logging_level FATAL" % (java, jar, gatk)
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            # CountReads
            #"-T CountReads -I /home/rob/spiral -R /home/rob/gatk-protected/public/gatk-engine/target/test-classes/exampleFASTA.fasta"

            # PrintReads
            self.check("SPEC_REF=%s %s -T PrintReads -I %s -R %s > %s.sam" % (ref, cmdline, golden_spec, golden_fasta, test_pat))
            self.check_results("%s.sam" % test_pat, golden_sam)

            # no SPEC_REF specified
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            cmd = "%s -T PrintReads -I %s -R %s > /dev/null 2>&1" % (cmdline, golden_spec, golden_fasta)
            self.assertEqual(1, subprocess.call(cmd, shell=True), "GATK should exit 1 with no SPEC_REF defined")

            # PrintReads with loci
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            self.check("SPEC_REF=%s %s -T PrintReads -I %s -R %s -L chr1:200 -L chr1:94610-94800 -o %s.bam" % (ref, cmdline, golden_spec, golden_fasta, test_pat))
            self.check_results("%s.bam" % test_pat, golden_loci_bam)

            # ClipReads: no CAT
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            self.check("SPEC_REF=%s %s -T ClipReads -I %s -R %s -X CAT -o %s.bam" % (ref, cmdline, golden_spec, golden_fasta, test_pat))
            self.check_results("%s.bam" % test_pat, golden_clipped_bam)

            # Encryption
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            test_spec = "%s/encrypted.spec" % self.data_dir
            self.check("SPEC_KEY=foobar bam2spec --in %s --out %s --ref %s --write-bai -f --no-match-reference >/dev/null" % (golden_bam, test_spec, ref))
            self.assertEqual(1, subprocess.call("SPEC_REF=%s %s -T PrintReads -I %s -R %s > %s.sam 2>/dev/null" % (ref, cmdline, test_spec, golden_fasta, test_pat), shell=True), "%s appears to be unencrypted" % test_spec)

            self.check("SPEC_KEY=foobar SPEC_REF=%s %s -T PrintReads -I %s -R %s > %s.sam" % (ref, cmdline, test_spec, golden_fasta, test_pat))
            self.check_results("%s.sam" % test_pat, golden_sam)


    @unittest.skipIf(not os.path.isfile('/opt/spiral/lib/spec-1.120.jar'), 'specstream not found')
    def test_zz_gatk_write(self):
        """
            spec write support for GATK
        """
        #self.cleanup = False

        ref = '/share/flat_refs/hg19.flat'

        golden_sam = '%s/spec/gatk/example_reads.sam' % GOLDEN_DIR
        golden_bam = '%s/spec/gatk/example_reads.bam' % GOLDEN_DIR
        golden_fasta = '%s/spec/gatk/exampleFASTA.fasta' % GOLDEN_DIR

        java_home = "/usr/lib/jvm/java-1.7.0-openjdk-amd64"
        java = "JAVA_HOME=%s %s/bin/java" % (java_home, java_home)
        gatk = "/share/gatk_resources/gatk-spiral-3.3-0.jar"

        test = 0
        for jar in ["/opt/spiral/lib/spec-1.120.jar", "/opt/spiral/lib/spec-1.130.jar"]:

            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            preamble = "SPEC_REF=%s %s -cp %s:%s org.broadinstitute.gatk.engine.CommandLineGATK" % (ref, java, jar, gatk)
            cmdline = "%s -T LeftAlignIndels -I %s -R %s -o %s.spec --logging_level FATAL" % (preamble, golden_bam, golden_fasta, test_pat)
            self.check(cmdline)
            self.check("spec2bam --in %s.spec --out %s.bam --ref %s" % (test_pat, test_pat, ref))
            self.check("samtools view -h %s.bam > %s.sam" % (test_pat, test_pat))
            self.check_results("%s.sam" % test_pat, golden_sam)

            # encryption
            test_pat = "%s/test%d" % (self.data_dir, test)
            test = test + 1

            preamble = "SPEC_KEY=foo SPEC_REF=%s %s -cp %s:%s org.broadinstitute.gatk.engine.CommandLineGATK" % (ref, java, jar, gatk)
            cmdline = "%s -T LeftAlignIndels -I %s -R %s -o %s.spec --logging_level FATAL" % (preamble, golden_bam, golden_fasta, test_pat)
            self.check(cmdline)

            # fail, wrong key
            self.assertEqual(2, subprocess.call("spec2bam --in %s.spec --out %s.bam --ref %s > /dev/null 2>&1" % (test_pat, test_pat, ref), shell=True), "%s.spec does not appear to be encrypted." % test_pat)

            self.check("SPEC_KEY=foo spec2bam --in %s.spec --out %s.bam --ref %s" % (test_pat, test_pat, ref))
            self.check("samtools view -h %s.bam > %s.sam" % (test_pat, test_pat))
            self.check_results("%s.sam" % test_pat, golden_sam)

    def test_spec_reference(self):
        """
            Test spec reference features
        """

        good_fasta = '{0}/spec/test1.fasta'.format(GOLDEN_DIR)
        bad_contig_fasta = '{0}/spec/test2.fasta'.format(GOLDEN_DIR)
        bad_md5_fasta = '{0}/spec/test3.fasta'.format(GOLDEN_DIR)
        bad_ids_fasta = '{0}/spec/broken.fasta'.format(GOLDEN_DIR)

        good_ref = '{0}/good_ref.flat'.format(self.data_dir)
        bad_contig_ref = '{0}/bad_contig_ref.flat'.format(self.data_dir)
        bad_md5_ref = '{0}/bad_md5_ref.flat'.format(self.data_dir)

        test_bam = '{0}/spec/test.bam'.format(GOLDEN_DIR)
        test_spec = '{0}/test.spec'.format(self.data_dir)
        test_local_bam = '{0}/test.bam'.format(self.data_dir)

        test_bam = '{0}/spec/test.bam'.format(GOLDEN_DIR)

        self.check("fasta2ref {0} {1}".format(good_fasta, good_ref))
        self.check("fasta2ref {0} {1}".format(bad_contig_fasta, bad_contig_ref))
        self.check("fasta2ref {0} {1}".format(bad_md5_fasta, bad_md5_ref))

        self.check("bam2spec --in {0} --ref {1} --out {2}  --no-match-reference > /dev/null 2>&1".format(test_bam, good_ref, test_spec))
        spec2bam_cmd = "spec2bam --in {0} --ref {1} --out {2} > /dev/null 2>&1".format(test_spec, good_ref, test_local_bam)
        self.check(spec2bam_cmd)
        os.remove(test_local_bam)
        spec2bam_cmd = "spec2bam --in {0} --ref {1} --out {2} > /dev/null 2>&1".format(test_spec, bad_contig_ref, test_local_bam)
        self.assertEqual(2, subprocess.call(spec2bam_cmd, shell=True), "spec2bam should fail on a reference with wrong contigs, but didn't")
        spec2bam_cmd = "spec2bam --in {0} --ref {1} --out {2} > /dev/null 2>&1".format(test_spec, bad_md5_ref, test_local_bam)
        self.assertEqual(2, subprocess.call(spec2bam_cmd, shell=True), "spec2bam should fail on a reference with different contigs, but didn't")

        cmd = "fasta2ref {0} invalid.ref > /dev/null 2>&1".format(bad_ids_fasta)
        self.assertEqual(2, subprocess.call(cmd, shell=True), "fasta2ref should fail on a reference with duplicate ids, but didn't")

        cmd = "bam2spec --in {0} --ref {1} --out {2}/none.spec -f > /dev/null 2>&1".format(test_bam, bad_contig_ref, self.data_dir)
        self.assertEqual(1, subprocess.call(cmd, shell=True), "bam2spec should fail on a reference with wrong contigs, but didn't")

        cmd = "bam2spec --in {0} --ref {1} --out {2}/none.spec -f --no-match-reference > /dev/null 2>&1".format(test_bam, bad_contig_ref, self.data_dir)
        self.check(cmd)

    # def test_spec2fasta(self):
    #     """
    #         Test spec to fasta conversion
    #     """
    #     golden_fasta = '{0}/spec/spec2fasta.fasta'.format(GOLDEN_DIR)
    #     golden_bam = '{0}/spec/test.bam'.format(GOLDEN_DIR)

    #     ref = '{0}/spec2fasta.ref'.format(self.data_dir)
    #     spec = '{0}/test.spec'.format(self.data_dir)
    #     fasta = '{0}/converted.fasta'.format(self.data_dir)

    #     self.check("fasta2ref --fasta-file {0} --ref-file {1}".format(golden_fasta, ref))
    #     self.check("bam2spec --in {0} --ref {1} --out {2} --embed-ref".format(golden_bam, ref, spec))

    #     self.check("spec2fasta {0} {1}".format(spec, fasta))
    #     self.check_results(golden_fasta, fasta)

    #     # Shouldn't work on a spec with no embedded ref
    #     self.check("bam2spec --in {0} --ref {1} --out {2} --force".format(golden_bam, ref, spec))
    #     self.assertEqual(2, subprocess.call("spec2fasta {0} {1} --force 2>/dev/null".format(spec, fasta), shell=True), "spec2fasta should exit non-zero when run on a spec with no embedded reference")

    def test_spec_full_roundtrip(self):
        """
            Test the entire spec round-trip
        """
        golden_fasta = '{0}/spec/spec2fasta.fasta'.format(GOLDEN_DIR)
        golden_bam = '{0}/spec/test.bam'.format(GOLDEN_DIR)

        ref = '{0}/spec2fasta.ref'.format(self.data_dir)
        spec = '{0}/test.spec'.format(self.data_dir)
        fasta = '{0}/converted.fasta'.format(self.data_dir)
        bam = '{0}/roundtrip.bam'.format(self.data_dir)
        decrypted_bam = '{0}/roundtrip.decrypted.bam'.format(self.data_dir)
        sam = '{0}/roundtrip.sam'.format(self.data_dir)

        self.check("fasta2ref {0} {1}".format(golden_fasta, ref))
        self.check("bam2spec --in {0} --ref {1} --out {2} --embed-ref --no-match-reference".format(golden_bam, ref, spec))

        self.check("bam2spec --in {0} --ref {1} --out {2}-fasta --embed-ref --no-match-reference".format(golden_bam, golden_fasta, spec))
        self.assertTrue(os.path.getsize(spec) == os.path.getsize("{0}-fasta".format(spec)), "Embedded fasta vs. ref are not the same size")

        # with crypto
        self.check("bam2spec --in {0} --ref {1} --out {2}.encrypted --embed-ref --key foo --no-match-reference".format(golden_bam, ref, spec))

        # use the embedded ref
        self.check("spec2bam --in {0} --out {1}".format(spec, bam))
        # use the external ref
        self.check("spec2bam --in {0}.encrypted --out {1} --ref {2} --key foo".format(spec, decrypted_bam, ref))

        self.check("samtools view {0} > {1}".format(bam, sam))
        self.check("samtools view {0} > {1}.decrypted".format(golden_bam, sam))
        self.check("samtools view {0} > {1}.golden".format(golden_bam, sam))

        self.check_results(sam, "{0}.golden".format(sam))
        self.check_results(sam, "{0}.decrypted".format(sam))

        self.check("spec2fasta {0} {1}".format(spec, fasta))
        self.check_results(golden_fasta, fasta)

    def test_spec_versions(self):
        """
            Refuse to work with incorrect versions
        """
        bad_major = '{0}/spec/test-99.1.0.spec'.format(GOLDEN_DIR)
        bad_minor = '{0}/spec/test-1.99.0.spec'.format(GOLDEN_DIR)

        ref = '/share/flat_refs/hg19.flat'

        cmd = "spec2bam --in {0} --ref {1} --out bad.bam --force >/dev/null 2>&1".format(bad_major, ref)
        self.assertEqual(2, subprocess.call(cmd, shell=True), "spec2bam should refuse to open files with the wrong major version")

        cmd = "spec2bam --in {0} --ref {1} --out bad.bam --force >/dev/null 2>&1".format(bad_minor, ref)
        self.assertEqual(2, subprocess.call(cmd, shell=True), "spec2bam should refuse to open files with the wrong minor version")

    def test_contig_remapping(self):
        """
            Contig remapping
        """
        bam = "{0}/spec/test.bam".format(GOLDEN_DIR)
        mapfile = "{0}/spec/remap.hg19_to_GRCh37.txt".format(GOLDEN_DIR)
        ref = "/share/flat_refs/homo_sapiens_GRCh37.flat"
        spec = "{0}/remapped.spec".format(self.data_dir)
        remapped = "{0}/remapped.bam".format(self.data_dir)
        reremapped = "{0}/reremapped.bam".format(self.data_dir)

        cmd = "bam2spec --in {0} --out {1} --ref {2} -f --remap {3}".format(bam, spec, ref, mapfile)
        self.check(cmd)
        cmd = "spec2bam --in {0} --out {1} --ref {2} -f --remap {3}".format(spec, remapped, ref, mapfile)
        self.check(cmd)
        cmd = "samtools view {0} |cut -f 3|sort -u".format(remapped)
        self.assertEqual("chr1", subprocess.check_output(cmd, shell=True).rstrip())

        cmd = "spec2bam --in {0} --out {1} --ref {2} -f".format(spec, reremapped, ref)
        self.check(cmd)
        cmd = "samtools view {0} |cut -f 3|sort -u".format(reremapped)
        self.assertEqual("1", subprocess.check_output(cmd, shell=True).rstrip())

    @unittest.skipUnless('PERFTEST' in os.environ, 'set PERFTEST=1 to enable')
    def test_spec_encode_perf(self):
        """
            encode perf testing
        """
        # self.cleanup = False
        golden_bam = '/share/bams/X10-sample.bam'
        ref = '/share/flat_refs/hg19.flat'
        spec_file = '%s/out.spec' % self.data_dir

        # Basic test
        print("", file=sys.stderr)
        for threads in [32, 8, 1]:
            start_time = time()
            self.check("bam2spec --in %s --ref %s --out %s --threads %d --force" % (golden_bam, ref, spec_file, threads))
            print("threads: %d elapsed time: %0.2f" % (threads, (time() - start_time)), file=sys.stderr)

    @unittest.skipUnless('PERFTEST' in os.environ, 'set PERFTEST=1 to enable')
    def test_spec_decode_perf(self):
        """
            decode perf testing
        """
        # self.cleanup = False
        golden_bam = '/share/bams/X10-sample.bam'
        ref = '/share/flat_refs/hg19.flat'
        spec_file = '%s/out.spec' % self.data_dir
        bam_file = '%s/out.bam' % self.data_dir

        self.check("bam2spec --in %s --ref %s --out %s --force" % (golden_bam, ref, spec_file))

        # Basic test
        print("", file=sys.stderr)
        for threads in [32, 8, 1]:
            start_time = time()
            self.check("spec2bamd --in %s --ref %s --out %s --threads %d --force" % (spec_file, ref, bam_file, threads))
            print("threads: %d elapsed time: %0.2f" % (threads, (time() - start_time)), file=sys.stderr)

if __name__ == '__main__':
    unittest.main(verbosity=2)
