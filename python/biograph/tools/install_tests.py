#!/usr/bin/env python3
"""
Tests to be run on an installation inside of lambdaToyData to ensure that
the installation was successful.
Checks that programs run correctly and that the SDK is available and that
query times are of the expected speed.
"""
# pylint: disable=invalid-name

from __future__ import print_function
import os
import sys
import time
import random
import unittest
import datetime
import argparse

try:
    import biograph
except ImportError as e:
    print("Error! BioGraph is not installed. Exiting")
    exit(1)

from biograph.utils import cmd_exe
import biograph.tools.log as log

def parse_args(args):
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog="discover", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--model", type=str, required=True,
                        help="BioGraph classifier model for qual_classifier")
    args = parser.parse_args(args)
    return args

def log_test(test_func):
    """
    Decorator to allow easy entering/exiting of tests to be written to the log
    """
    def wrapper(*args):
        """
        Wrapper
        """
        log.info("Checking %s", test_func.__name__)
        test_func(*args)
        log.info("Completed %s", test_func.__name__)
    return wrapper


class BioGraphInstallTest(unittest.TestCase):

    """
    Tests to evaluate if installing biograph was successful
    """

    ML_MODEL_FN = None

    def setUp(self):
        """
        Ensures we're in the tarball
        """
        # Lookups for the sample: read_count, base_count, query's read count
        self.query_seq = "CATCATCAT"
        self.rm_cnts = {"proband": [48956, 7235142, 122],
                        "father": [48930, 7228335, 124],
                        "mother": [48964, 7235912, 151]}

    def create_check(self, ret, out, err, tdel):
        """
        Do the work of checking per individual
        """
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertEqual(ret, 0, "BioGraph sample creation returned non-zero")
        self.assertTrue(tdel < datetime.timedelta(minutes=4),
                        "Creation time took longer than 4 minutes (%s)" % str(tdel))

    def sdk_file_structure_work(self, bgname):
        """
        Using the sdk to ensure that the BioGraphs were all created correctly
        by looking at the metadata and checking on the file structure
        """
        log.info("Opening %s", bgname)
        my_bg = biograph.BioGraph(bgname)
        log.info("BioGraph Id: %s", my_bg.metadata.biograph_id)
        log.info("Accession Id: %s", my_bg.metadata.accession_id)
        self.assertTrue(os.path.exists(os.path.join(bgname, 'seqset')))
        for sample in my_bg.metadata.samples:
            my_rm = my_bg.open_readmap(sample)
            readmap = os.path.join(bgname, 'coverage',
                                   my_bg.metadata.samples[sample] + '.readmap')
            self.assertTrue(os.path.exists(readmap), "Readmap data for sample %s does not exist" % (sample))
            self.assertEqual(my_rm.get_read_count(), self.rm_cnts[sample][0],
                             "Unexpected number of reads in %s BioGraph (%d)" % (sample, my_rm.get_read_count()))
            self.assertEqual(my_rm.get_num_bases(), self.rm_cnts[sample][1],
                             "Unexpected number of bases in %s BioGraph (%d)" % (sample, my_rm.get_num_bases()))
            read_cnt = len(list(my_rm.get_reads_containing(my_bg.seqset.find(self.query_seq))))
            if sample != 'father':  # Edge case for merge
                self.assertEqual(read_cnt, self.rm_cnts[sample][2],
                                 "Unexpected number of reads containing '%s' in %s BioGraph (%s)"
                                 % (self.query_seq, sample, read_cnt))

    @log_test
    def test000_ensure_directory(self):
        """
        First, ensure that all the files are in the places we expect them to be
        """
        self.assertTrue(os.path.exists("../references/reference_lambda.fasta"), "Reference file is missing")
        self.assertTrue(os.path.exists("../reads/proband_1.fastq.gz"), "Proband read1 fastq is missing")
        self.assertTrue(os.path.exists("../reads/proband_2.fastq.gz"), "Proband read2 fastq is missing")
        self.assertTrue(os.path.exists("../reads/father_reads.bam"), "Father reads bam is missing")
        self.assertTrue(os.path.exists("../reads/mother_reads.bam"), "Mother reads bam is missing")
        self.assertTrue(os.path.exists("../pedigree.ped"), "pedigree is missing")
        self.assertTrue(os.path.exists("../variants/family.vcf.gz"), "Family vcf.gz is missing")
        self.assertTrue(os.path.exists("../variants/family.vcf.gz.tbi"), "Family vcf.gz.tbi is missing")

    @log_test
    def test010_reference(self):
        """
        Creates the lambda reference
        """
        ret, out, err, tdel = cmd_exe(("biograph reference "
                                       "--in ../references/reference_lambda.fasta "
                                       "--refdir ref_lambda"))
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertTrue(os.path.exists('ref_lambda'), 'Reference folder does not exist')
        self.assertTrue(tdel < datetime.timedelta(minutes=2),
                        "Reference creation time took longer than 2 minutes (%s)" % str(tdel))

    @log_test
    def test015_reference_sdk(self):
        """
        Tests the Reference SDK and validity
        """
        ref = biograph.Reference("ref_lambda")
        match = ref.find("AGTACAT")
        self.assertEqual(match.matches, 1, "Unexpected number of matches over test query")
        match = match.get_match(0)
        self.assertEqual(match.chromosome, 'lambda', "Unexpected chromosome found (%s)" % match.chromosome)
        self.assertEqual(match.start, 29071, "Unexpected reference match start found (%d)" % match.start)
        self.assertEqual(match.end, 29078, "Unexpected reference match start found (%d)" % match.end)

    @log_test
    def test020_create_proband(self):
        """
        Creates the proband biograph and ensures biograph was run correctly
        """
        log.info("Creating Proband BioGraph")
        self.create_check(*cmd_exe(("biograph create "
                                    "--in ../reads/proband_1.fastq.gz "
                                    "--pair ../reads/proband_2.fastq.gz "
                                    "--ref ref_lambda/ "
                                    "--out proband_lambda.bg "
                                    "--id proband")))

    @log_test
    def test030_sdk_file_structure_proband(self):
        """
        Checks results are there
        """
        log.info("Checking proband.bg file structure")
        self.sdk_file_structure_work("proband_lambda.bg")
        my_bg = biograph.BioGraph("proband_lambda.bg")
        ssz = my_bg.seqset.size()
        self.assertEqual(ssz, 97340, "Unexpected Seqset size (%d)" % (ssz))

    @log_test
    def test040_create_father(self):
        """
        Creates the father biograph and ensures biograph was run correctly
        """
        log.info("Creating Father BioGraph")
        self.create_check(*cmd_exe(("biograph create "
                                    "--in ../reads/father_reads.bam "
                                    "--ref ref_lambda/ "
                                    "--out father_lambda.bg "
                                    "--id father")))

    @log_test
    def test050_sdk_file_structure_father(self):
        """
        Checks results are there
        """
        log.info("Checking father_lambda.bg file structure")
        self.sdk_file_structure_work("father_lambda.bg")
        my_bg = biograph.BioGraph("father_lambda.bg")
        ssz = my_bg.seqset.size()
        self.assertEqual(ssz, 98006, "Unexpected Seqset size (%d)" % (ssz))

    @log_test
    def test060_create_mother(self):
        """
        Creates the mother biograph and ensures biograph was run correctly
        """
        log.info("Creating Mother BioGraph")
        self.create_check(*cmd_exe(("biograph create "
                                    "--in ../reads/mother_reads.bam "
                                    "--ref ref_lambda/ "
                                    "--out mother_lambda.bg "
                                    "--id mother")))

    @log_test
    def test070_sdk_file_structure_mother(self):
        """
        Checks results are there
        """
        log.info("Checking mother.bg file structure")
        self.sdk_file_structure_work("mother_lambda.bg")
        my_bg = biograph.BioGraph("mother_lambda.bg")
        ssz = my_bg.seqset.size()
        self.assertEqual(ssz, 96697, "Unexpected Seqset size (%d)" % (ssz))

    @unittest.skip("New biograph wrapper script doesn't support merging")
    @log_test
    def test080_merge(self):
        """
        Test that merging is functioning correctly
        """
        log.info("Merging Three BioGraphs")
        ret, out, err, tdel = cmd_exe(("biograph merge "
                                       "--in proband_lambda.bg/ "
                                       "--in father_lambda.bg/ "
                                       "--in mother_lambda.bg/ "
                                       "--out family_lambda.bg "
                                       "--id family"))
        self.assertEqual(ret, 0, "BioGraph merging returned non-zero")
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertTrue(tdel < datetime.timedelta(minutes=1),
                        "Creation time took longer than 1 minute (%s)" % str(tdel))

    # @log_test
    # def test090_sdk_file_structure_merge(self):
    #     """
    #     Test the merged file structure and that all the samples are correct
    #     """
    #     log.info("Checking family_lambda.bg file structure")
    #     self.sdk_file_structure_work("family_lambda.bg")
    #     my_bg = biograph.BioGraph("family_lambda.bg")
    #     ssz = my_bg.seqset.size()
    #     self.assertEqual(ssz, 103996, "Unexpected Seqset size (%d)" % (ssz))
    #     cnt = len(my_bg.metadata.samples.keys())
    #     self.assertEqual(cnt, 3, "Unexpected number of Readmaps in merged BioGraph (%d)" % (cnt))

    # @log_test
    # def test100_sdk_query_time(self):
    #     """
    #     Using the sdk to ensure that the BioGraph queries perform as well as we expect them to
    #     """
    #     my_bg = biograph.BioGraph("proband_lambda.bg")
    #     ref = biograph.Reference("ref_lambda")
    #     start_time = time.time()
    #     _ = biograph.find_region_variants(my_bg, ref, "lambda", 0, 48503)
    #     end_time = time.time()
    #     exec_time = datetime.timedelta(seconds=int(end_time - start_time))
    #     self.assertTrue(exec_time < datetime.timedelta(seconds=1),
    #                     "find_region_variants took more than 1 seconds (%s)" % str(exec_time))
    #     start_time = time.time()
    #     _ = my_bg.get_reads_containing("CAGCAGCACAG", 2, 1000)
    #     end_time = time.time()
    #     exec_time = datetime.timedelta(seconds=int(end_time - start_time))
    #     self.assertTrue(exec_time < datetime.timedelta(seconds=1),
    #                     "find_region_variants took more than 1 seconds (%s)" % str(exec_time))

    @log_test
    def test110_bgtools_works(self):
        """
        See if bgtools is installed correctly and we can use coverage to annotate the vcf
        """
        ret, out, err, tdel = cmd_exe("biograph -h")
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertEqual(ret, 0, "There was a problem running `biograph`")
        start_time = time.time()
        ret, out, err, tdel = cmd_exe(("biograph coverage -b proband_lambda.bg "
                                       "-v ../variants/family.vcf.gz "
                                       "-r ref_lambda/ "
                                       "-s proband "
                                       "-o proband_pcmp.vcf"))
        end_time = time.time()
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertEqual(ret, 0, "There was a problem running `biograph coverage`")
        exec_time = datetime.timedelta(seconds=int(end_time - start_time))
        self.assertTrue(exec_time < datetime.timedelta(seconds=5),
                        "biograph coverage took more than 5 seconds (%s)" % str(exec_time))

    @log_test
    def test120_full_pipeline(self):
        """
        Make sure the full pipeline can be run.
        """
        start_time = time.time()
        ret, out, err, tdel = cmd_exe(
            "biograph full_pipeline"
            " --biograph full_pipeline_father.bg"
            " --reference ref_lambda"
            " --reads ../reads/father_reads.bam"
            " --model " + BioGraphInstallTest.ML_MODEL_FN)
        end_time = time.time()
        log.info("RETCODE: %d", ret)
        log.info("CMDTIME: %s", str(tdel))
        log.info("STDOUT: %s", out)
        log.info("STDERR: %s", err)
        self.assertEqual(ret, 0, "There was a problem running `biograph full_pipeline`")
        exec_time = datetime.timedelta(seconds=int(end_time - start_time))
        self.assertTrue(exec_time < datetime.timedelta(seconds=60),
                        "biograph full pipeline took more than 60 seconds (%s)" % str(exec_time))



class cd:

    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
        self.savedPath = None

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class LogFileStdout:

    """ Write to stdout and a file"""

    def __init__(self, fn):
        """ keep these props """
        self.name = fn
        self.file_handler = open(fn, 'w')

    def write(self, *args):
        """ Write to both """
        sys.stdout.write(*args)
        self.file_handler.write(*args)

    def flush(self):
        """ Flush both """
        self.file_handler.flush()
        sys.stdout.flush()


def main(args):
    """
    Runs the tests for validating the installation
    """
    args = parse_args(args)
    ret, unused, unused, unused = cmd_exe("bgbinary")
    if ret == 127:
        print("Error! The BioGraph executable is not in path")
        exit(1)
        print(unused)
    start_time = time.time()
    out_dir = ("test_%8x" % random.randrange(16 ** 8)).replace(' ', '_').replace('\n', '_')
    os.mkdir(out_dir)
    print("Writing test results to %s" % os.path.abspath(out_dir))
    tprog = None
    with cd(out_dir):
        log_out = LogFileStdout('log.txt')
        log.setup_logging(stream=log_out)
        print("See %s/%s for details on the run" % (out_dir, log_out.name))
        unittest.TestLoader.sortTestMethodsUsing = None
        BioGraphInstallTest.ML_MODEL_FN = args.model
        suite = unittest.TestLoader().loadTestsFromTestCase(BioGraphInstallTest)
        tprog = unittest.TextTestRunner(log_out, verbosity=2, failfast=False).run(suite)
    end_time = time.time()
    exec_time = datetime.timedelta(seconds=int(end_time - start_time))
    log.info("Tests Completed in %s", str(exec_time))
    if tprog.wasSuccessful():
        print("BioGraph has been successfully installed!")
    else:
        print("There were issues validating the BioGraph installation.")
        print("Check %s/log.txt and contact Spiral for assistance." % os.path.abspath(out_dir))
        exit(1)

if __name__ == '__main__':
    main(sys.argv[1:])
