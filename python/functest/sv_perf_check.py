"""
sv_performance_test.py: Medium tests to ensure variant discovery hasn't changed too much
"""

import os
import json
import unittest
import sys
import subprocess

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)

SEARCH_PATH = ["/scratch"]
if "HOME" in os.environ:
    SEARCH_PATH.append(os.environ["HOME"] + "/datasets")
if "USER" in os.environ:
    SEARCH_PATH.append("/home/" + os.environ["USER"] + "/datasets")
SEARCH_PATH = SEARCH_PATH + ["/share/datasets", "/share/datasets/HG002"]

REFMAP_DIR = "/scratch"

def get_dataset_path(filename):
    """Searches for filename in SEARCH_PATH"""
    for p in SEARCH_PATH:
        candidate = os.path.join(p, filename)
        if os.path.exists(candidate):
            print("Found %s in %s" % (filename, candidate))
            sys.stdout.flush()
            return candidate
    raise FileNotFoundError("Could not find %s in search path" % filename)

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()


class SVBenchmark(unittest.TestCase):
    """ unittest Test definitions follow """
    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def test_00_create_svs(self): # pylint: disable=too-many-locals
        """
        This is the full test
        """
        os.environ['PATH'] = os.environ['PATH'] + ":/share/software/bin"

        bgname = "HG002-NA24385-50x.bg"
        bedfile = get_dataset_path("HG002/giab_v0.6/HG002_SVs_Tier1_v0.6.chr20.bed")
        # bedfile = get_dataset_path("mini-chr20.bed")
        truth_v = get_dataset_path("HG002/giab_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz")
        bg_file = get_dataset_path("HG002/" + bgname)
        ref_seq = "/reference/hs37d5/"
        truth_baseline = "python/functest/sv_perf_baseline.txt"
        out_var = os.path.join(self.data_dir, "var.vcf")

        refmap_arg = ""
        if os.access(REFMAP_DIR, os.W_OK):
            refmap_arg = " --ref-map " + REFMAP_DIR + "/" + bgname + ".refmap"
        bg_cmd = "bgbinary discovery --in {0} --bed {1} --ref {2} --out {3} {4}".format(bg_file, bedfile, ref_seq, out_var, refmap_arg)
        self.assertEqual(0, subprocess.call(bg_cmd, shell=True))

        zip_out = out_var + ".gz"
        sort_cmd = "vcf-sort {0} | bgzip > {1} && tabix {1}".format(out_var, zip_out)
        self.assertEqual(0, subprocess.call(sort_cmd, shell=True))

        tru_out = os.path.join(self.data_dir, "truvari_sv")
        tru_cmd = "/usr/bin/env python3.6 external/*pypi__Truvari_*/Truvari-*.data/scripts/truvari -b {0} -c {1} -f {2}/source.fasta -o {3} --passonly --includebed {4}".format(
            truth_v, zip_out, ref_seq, tru_out, bedfile)
        self.assertEqual(0, subprocess.call(tru_cmd, shell=True))

        with open(truth_baseline, 'r') as f:
            base = json.load(f)
        with open(os.path.join(tru_out, "summary.txt"), 'r') as f:
            comp = json.load(f)
        p_diff = comp["precision"] - base["precision"]
        r_diff = comp["recall"] - base["recall"]

        fn_sample = subprocess.check_output(f"egrep -v '^#' {tru_out}/fn.vcf | shuf -n 10", shell=True).decode()
        err = f"""
Precision went from {base['precision']*100:.2f}% to {comp['precision']*100:.2f}%, a {p_diff*100:.2f}% difference
Recall went from {base['recall']*100:.2f}% to {comp['recall']*100:.2f}%, a {r_diff*100:.2f}% difference

Random sample of false negatives:
{fn_sample}
"""

        self.assertLess(abs(p_diff), .03, err)
        self.assertLess(abs(r_diff), .03, err)

        # Still print the stats for informational purposes, even if we succeeded"
        print(err)

if __name__ == '__main__':
    unittest.main(verbosity=2)
