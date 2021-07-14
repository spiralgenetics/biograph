"""
sdk.py: Test the python SDK and ipython notebooks.
"""
import os
import glob
import queue
import unittest
import subprocess

import tabix

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

import biograph
from biograph.tools.coverage import MICROCONTIGS, get_regions

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()


class SDKTestCases(unittest.TestCase):

    """    unittest Test definitions follow """

    # keep pylint happy
    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit code is correct """
        actual = subprocess.call(cmd # + ">/dev/null 2>&1"
                                 , shell=True)
        self.assertEqual(code, actual, f"{cmd} exited code {actual}")

    def check_results(self, file1, file2):
        """ sha1 two files, throw if they don't match """
        self.assertEqual(sha1_file(test=self, file_name=file1),
                         sha1_file(file2), "%s and %s do not match" % (file1, file2))

    # @unittest.skip('nope')
    def test_basic_sdk(self):
        """
            Basic SDK test
        """
        my_bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/proband_lambda.bg")

        self.assertEqual(my_bg.metadata.accession_id, 'proband')
        self.assertEqual(my_bg.metadata.samples['proband'], '47b43dbb0306fb5da00bca9c33b2fb3de9db7bf4')
        self.assertEqual(str(my_bg.seqset.find('ACGT')), '<SeqsetEntry 11049-11336: ACGT>')
        self.assertEqual(my_bg.seqset.size(), 97340)

        r = my_bg.open_readmap()
        self.assertEqual(r.get_num_bases(), 7235142)

        stats = r.get_pair_stats()
        self.assertEqual(stats.paired_reads, 39974)
        self.assertEqual(stats.unpaired_reads, 8982)


    # @unittest.skip('nope')
    def test_basic_sdk_reference(self):
        """
            Basic reference test
        """
        grch37 = biograph.Reference("/reference/human_g1k_v37/")

        # 16 contigs in chr 7
        self.assertEqual(len(grch37.find_ranges("7", 1, 1000000000)), 16)

        # 11 contigs in chr 21
        self.assertEqual(len(grch37.find_ranges("21", 1, 1000000000)), 11)

        rc = grch37.find('ACGT')
        self.assertEqual(rc.matches, 2153084)

    # @unittest.skip('nope')
    def test_sequence_splice(self):
        """
        Testing the splicing/indexing and string comparisons for sequence
        """
        string = "ATCGAATACATAA"

        j = biograph.Sequence(string)
        # Some random splicing and access
        self.assertEqual(string[1], j[1])
        self.assertEqual(string[1:5], j[1:5])
        # self.assertEqual(string[2:7:2], j[2:7:2])
        # self.assertEqual(string[::-1], j[::-1])

    # @unittest.skip('nope')
    def test_bgtools_basic(self):
        """ bgtools script tests """
        bgtools = "./python/functest/biograph_main "

        # exit 1
        self.check(bgtools, code=1)

        # exit 0
        self.check(f"{bgtools} -h")

        # version
        cmd = f"{bgtools} version"
        self.assertEqual(subprocess.check_output(cmd, shell=True).rstrip().decode(), f"biograph version {biograph.version()}")

        # Every command we ship should have a help
        for cmd in (
                "full_pipeline",
                "reference",
                "create",
                "discovery",
                "coverage",
                "qual_classifier",
                "vdb",
                "stats",
            ):
            self.check(f"{bgtools} {cmd} -h")

    # @unittest.skip('nope')
    def test_bgtools_stats(self):
        """ bgtools script tests """
        self.maxDiff = None
        bgtools = "./python/functest/biograph_main "

        cmd = f"{bgtools} stats -b golden/e_coli_merged.bg -r datasets/reference/e_coli_k12_ASM584v1/"
        # pylint: disable=bad-continuation
        self.assertEqual(subprocess.check_output(cmd, shell=True).decode(),
"""Sample:            e_coli_test
NumReads:          38,047
NumBases:          3,804,700
MaxReadLength:     100
MinReadLength:     100
NumPairedReads:    25,762
NumUnpairedReads:  12,285
NumPairedBases:    2,576,200
NumUnpairedBases:  1,228,500
MeanInsertSize:    504.73
MedianInsertSize:  502.00
SDInsertSize:      49.72
EstimatedCoverage: 0.80

Sample:            test_accession_id
NumReads:          8,444
NumBases:          288,464
MaxReadLength:     35
MinReadLength:     30
NumPairedReads:    0
NumUnpairedReads:  8,444
NumPairedBases:    0
NumUnpairedBases:  288,464
MeanInsertSize:    0.00
MedianInsertSize:  0.00
SDInsertSize:      0.00
EstimatedCoverage: 0.06

""")

    # @unittest.skip('nope')
    def test_bgtools_pcmp(self):
        """ bgtools script tests """
        os.environ['PATH'] = "/share/software/bin:" + os.environ['PATH']

        bgtools = "./python/functest/biograph_main "
        bgb = "./modules/biograph/bgbinary "
        dd = self.data_dir

        self.check(f"{bgb} create --in datasets/bams/e_coli/e_coli_test.bam --out {dd}/e_coli.bg --ref datasets/reference/e_coli_k12_ASM584v1/")
        self.check(f"{bgb} discovery --in {dd}/e_coli.bg --out {dd}/e_coli_variants.vcf --enable-pop-tracer=false --ref datasets/reference/e_coli_k12_ASM584v1/")
        self.assertEqual(subprocess.call(f"vcf-sort < {dd}/e_coli_variants.vcf | bgzip > {dd}/e_coli_variants.vcf.gz", shell=True), 0)
        self.check(f"tabix {dd}/e_coli_variants.vcf.gz")
        self.check(f"{bgtools} coverage -d {dd}/ecoli.cov.jl -b {dd}/e_coli.bg -v {dd}/e_coli_variants.vcf.gz -r datasets/reference/e_coli_k12_ASM584v1/ -o /dev/stdout | vcf-sort > {dd}/output.vcf")

        self.check(f"{bgb} discovery --in {dd}/e_coli.bg --out {dd}/e_coli_variants_ml.vcf --ref datasets/reference/e_coli_k12_ASM584v1/")
        self.assertEqual(subprocess.call(f"vcf-sort < {dd}/e_coli_variants_ml.vcf | bgzip > {dd}/e_coli_variants_ml.vcf.gz", shell=True), 0)
        self.check(f"tabix {dd}/e_coli_variants_ml.vcf.gz")
        self.check(f"{bgtools} coverage -d {dd}/e_coli_ml.cov.jl -b {dd}/e_coli.bg -v {dd}/e_coli_variants_ml.vcf.gz -r datasets/reference/e_coli_k12_ASM584v1/ -o /dev/stdout | vcf-sort > {dd}/output2.vcf")

        # Better test would be truvari to validate the actual variants, but an identical output line count is close enough for now.
        self.assertEqual(subprocess.check_output(f"wc -l < {dd}/output.vcf", shell=True).rstrip(), b'166', 0)

        # This will exit 1 if exceptions correctly propagate from child workers back up to the main process (DEV-511)
        self.check(f"{bgtools} coverage -d {dd}/ecoli.cov.jl -b {dd}/e_coli.bg -v golden/ftest/bad.vcf.gz -r datasets/reference/e_coli_k12_ASM584v1 -o {dd}/bad_out.vcf", 1)


    # @unittest.skip('nope')
    def test_bgtools_qual_classifier(self):
        """ bgtools script tests """
        os.environ['PATH'] = "/share/software/bin:" + os.environ['PATH']
        dd = self.data_dir

        bgtools = "./python/functest/biograph_main "
        toydata = "datasets/ml_toydata_lambda"
        sample_dir = "proband_jul17"
        golden_file = "proband_17feb2021"
        df = f"{toydata}/{sample_dir}/df.jl"
        grm = f"{toydata}/{sample_dir}/grm.jl"
        pcmp = f"{toydata}/{sample_dir}/pcmp.vcf.gz"
        sv_ao_model = "python/functest/biograph_model.ml"

        # Strip headers and compare variants for sv and ao
        self.check(f"{bgtools} qual_classifier --grm {grm} --dataframe {df} --vcf {pcmp} --model {sv_ao_model} --out {dd}/sv_ao.vcf")
        self.check(f"grep -v ^# {GOLDEN_DIR}/ml/{golden_file}.filter.vcf > {dd}/sv_ao.golden.vcf")
        self.check(f"vcf-sort < {dd}/sv_ao.vcf | grep -v ^# > {dd}/sv_ao.check.vcf")
        self.check_results(f"{dd}/sv_ao.golden.vcf", f"{dd}/sv_ao.check.vcf")


    # @unittest.skip('nope')
    def test_coverage_microcontigs(self):
        """ microcontig boundary generation tests """

        contigq = queue.Queue()
        workerq = queue.Queue()
        variants = set()

        ref = biograph.Reference('datasets/reference/e_coli_k12_ASM584v1/')

        # This VCF should cover every edge case (multiple variants on the same
        # position, deletions that span a boundary, begin, end, one before begin, one
        # after end, with and without a bed file, etc.)
        #
        # This tests counts on the vcf ID field being unique for each variant
        vcf_file = f'{GOLDEN_DIR}/microcontigs.vcf.gz'

        for bed_file in None, f'{GOLDEN_DIR}/microcontigs.bed':
            for clearance in range(2, 10):
                for contig_size in range(1, 20):
                    for region in get_regions(bed_file, ref):
                        contigq.put(region)
                    contigq.put(None)

                    # Find the microcontigs
                    while True:
                        item = contigq.get()
                        if item is None:
                            break
                        MICROCONTIGS.find(item, vcf_file, clearance, contig_size, workerq)

                    workerq.put(None)

                    # Add every variant
                    t = tabix.open(vcf_file)
                    for ctg in ref.scaffolds:
                        try:
                            for v in t.query(ctg, 0, int(ref.scaffold_lens[ctg])):
                                variants.add(v[2])
                        except tabix.TabixError:
                            pass

                    self.assertTrue(len(variants) > 0, f"No variants to process. Are you using the correct reference? (bed: {bed_file} clearance: {clearance} contig_size: {contig_size})")

                    # Remove each variant, one microcontig at a time
                    while True:
                        contig = workerq.get()
                        if contig is None:
                            break
                        # print(contig)
                        for v in t.query(contig[0], contig[1], contig[2]):
                            # This will raise if the variant was previously removed
                            variants.remove(v[2])

                    # There should be none left over
                    self.assertEqual(len(variants), 0, f"There were {len(variants)} unprocessed variants remaining.  (bed: {bed_file} clearance: {clearance} contig_size: {contig_size})")

    def setup_dirs(self):
        "Makes a bin directory with bgbinary, biograph, and truvari, and adds it to the path."
        bindir = os.path.join(self.data_dir, "bin")
        tmpdir = os.path.join(self.data_dir, "tmp")
        for dirname in bindir, tmpdir:
            os.mkdir(dirname)
        os.symlink(os.getcwd() + "/modules/biograph/bgbinary", bindir + "/bgbinary")
        truvari_dir = glob.glob(os.getcwd() + "/external/*pypi__Truvari_2*/Truvari-*.data/scripts")[0]
        os.symlink(os.getcwd() + "/python/functest/biograph_main", bindir + "/biograph")
        with open(bindir + "/truvari", "w") as f:
            f.write(f"""#!/bin/sh\nexec python3 {truvari_dir}/truvari "$@"\n""")
        os.chmod(bindir + "/truvari", 0o755)

        os.environ['PATH'] = (bindir + ":" +
                              os.environ['PATH'] + ":" +
                              "/share/software/bin")
        return tmpdir



    @unittest.skip('lambda toydata is broken')
    def test_bgtools_squareoff(self):
        """ bgtools script tests """
        self.setup_dirs()

        dd = self.data_dir

        toydata = "/share/ml_toydata_lambda"
        sample_dir = "proband_jul17"
        golden_file = "proband_16feb2021"
        pcmp = f"{toydata}/{sample_dir}/pcmp.vcf.gz"
        sv_ao_model = "/share/classifier/model_sv_ao_gt_16Feb2021.ml"
        bg = "datasets/lambdaToyData/benchmark/proband_lambda.bg"
        ref = "datasets/lambdaToyData/benchmark/ref_lambda/"

        self.check(f"biograph squareoff -b {bg} --variants {pcmp} --model {sv_ao_model} --reference {ref} -o {dd}/sq.vcf.gz")
        self.check(f"zcat {dd}/sq.vcf.gz | vcf-sort | grep -v ^# > {dd}/sq.check.vcf")
        self.check(f"grep -v ^# {GOLDEN_DIR}/ml/{golden_file}.filter.vcf > {dd}/sv_ao.golden.vcf")
        with open(f"{dd}/sv_ao.golden.vcf") as fha, open(f"{dd}/sq.check.vcf") as fhb:
            acnt = len(fha.readlines())
            bcnt = len(fhb.readlines())
            self.assertEqual(acnt, bcnt, f"Incorrect line counts {acnt} != {bcnt}")

#     @unittest.skip('runipy output is not deterministic')
#     def test_notebooks(self):
#         """
#             Test all jupyter (iPython) notebooks

#             This doesn't currently work because runipy doesn't produce deterministic output
#             when the script writes to STDOUT (the output fields are arbitrarily chunked and
#             impossible to compare without further manipulation.)

#             To make golden files for new notebooks, run this:

#               runipy new_notebook.ipynb $GOLDEN_DIR/jupyter/new_notebook.ipynb.out
#         """
#         ipydir = "{0}/python/jupyter/biograph/".format(SPIRAL_ROOT)
#         for f in os.listdir(ipydir):
#             if f.endswith(".ipynb"):
#                 self.check("runipy '{0}/{1}' '{2}/{1}.out'".format(ipydir, f, self.data_dir))
#                 # strip out "matplotlib.figure.Figure at <0xaddress>" from txt lines
#                 subprocess.call(
#                     "grep -v matplotlib.figure.Figure '{0}/{1}.out' > '{0}/{1}.out.clean'".format(self.data_dir, f), shell=True)
#                 self.check_results('{0}/jupyter/{1}.out'.format(
#                     GOLDEN_DIR, f), '{0}/{1}.out.clean'.format(self.data_dir, f))


if __name__ == '__main__':
    unittest.main(verbosity=2)
