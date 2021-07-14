"""
sv_performance_test.py: Medium tests to ensure variant discovery hasn't changed too much
"""

import os
import json
import unittest
import sys
import subprocess
from glob import glob

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

BG_VERSION = sys.argv[1]
sys.argv = sys.argv[1:]
BED = "python/functest/sv_perf_check_full.bed"
BAM = "/share/datasets/HG002/HG002.hs37d5.10x.bam"
MODEL = "/share/releases/biograph_model-" + BG_VERSION + ".ml"
OUTPUT_TYPES = ("discovery", "results")

# Use this biograph instead of running the full pipeline.
# useful for testing this test without it taking so long:
SKIP_FULL_PIPELINE_RUN = None

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

    def setup_dirs(self):
        "Makes a bin directory with bgbinary, biograph, and truvari, and adds it to the path."
        bindir = os.path.join(self.data_dir, "bin")
        tmpdir = os.path.join(self.data_dir, "tmp")
        for dirname in bindir, tmpdir:
            os.mkdir(dirname)
        os.symlink(os.getcwd() + "/modules/biograph/bgbinary", bindir + "/bgbinary")
        truvari_dir = glob(os.getcwd() + "/external/*pypi__Truvari_2*/Truvari-*.data/scripts")[0]
        os.symlink(os.getcwd() + "/python/functest/biograph_main", bindir + "/biograph")
        with open(bindir + "/truvari", "w") as f:
            f.write(f"""#!/bin/sh
exec python3 {truvari_dir}/truvari "$@"
""")
        os.chmod(bindir + "/truvari", 0o755)

        os.environ['PATH'] = (bindir + ":" +
                              os.environ['PATH'] + ":" +
                              "/share/software/bin")
        return tmpdir

    def test_00_create_svs(self): # pylint: disable=too-many-locals
        """
        This is the full test
        """
        tmpdir = self.setup_dirs()
        biograph_dir = os.path.join(self.data_dir, "sv_perf_check_full.bg")
        bg_cmd = f"""biograph full_pipeline
        -b {biograph_dir}
        -r /reference/hs37d5
        --reads {BAM}
        --model {MODEL}
        --create '--max-mem 60 --min-kmer-count=2'
        --discovery '--bed {BED}'
        --coverage '-R {BED}'
        --tmp {tmpdir}
        --keep all""".replace("\n", "")

        if SKIP_FULL_PIPELINE_RUN:
            bg_cmd = f"ln -s {SKIP_FULL_PIPELINE_RUN} {biograph_dir}"
        self.assertEqual(0, subprocess.call(bg_cmd, shell=True))

        truth_v = get_dataset_path("HG002/giab_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz")
        ref_seq = "/reference/hs37d5/"
        truth_baseline = "python/functest/sv_perf_full_baseline.txt"

        passed = True
        try:
            with open(truth_baseline, 'r') as f:
                base = json.load(f)
        except Exception as e: # pylint: disable=broad-except
            print(f"ERROR: No base data available: {e}")
            base = None
            passed = False

        truvaris = {}

        for output_type in OUTPUT_TYPES:
            tru_out = os.path.join(self.data_dir, "truvari_sv_" + output_type)
            out_var = f"{biograph_dir}/analysis/{output_type}.vcf.gz"
            tru_cmd = f"/usr/bin/env python3.6 external/*pypi__Truvari_1_1/Truvari-*.data/scripts/truvari -b {truth_v} -c {out_var} -f {ref_seq}/source.fasta -o {tru_out} --passonly --includebed {BED}"
            self.assertEqual(0, subprocess.call(tru_cmd, shell=True))
            with open(os.path.join(tru_out, "summary.txt"), 'r') as f:
                truvaris[output_type] = json.load(f)
        print(f"Truvari evaluations:\n{json.dumps(truvaris)}\n")


        for output_type in OUTPUT_TYPES:
            comp = truvaris[output_type]
            for stat_type in ("precision", "recall"):
                comp_val = comp[stat_type]
                if base:
                    base_val = base[output_type][stat_type]
                    diff = comp_val - base_val
                    if abs(diff) < 0.03:
                        status = "OK"
                    else:
                        status = "FAIL"
                        passed = False
                else:
                    passed = False
                    status = "UNKNOWN"
                    base_val = "UNKNOWN"
                    diff = 0
                print(f"{status}: {output_type} {stat_type} = {comp_val}, base = {base_val}, diff={diff*100.0:.2f}%")

        self.assertTrue(passed)

if __name__ == '__main__':
    unittest.main(verbosity=2)
