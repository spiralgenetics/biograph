#!/usr/bin/env python3
"""
vdb_test.py: Test the Variant DB. Run manually because these tests cost money.

To run under bazel, using the default "vdb" database and "spiral-vdb" bucket:
  bazel test //python/functest:vdb_test
"""
import unittest
import subprocess
import tempfile
import os
import logging

# from pathlib import Path
from time import time

from python.functest.utils.fileops import sha1_file

import biograph.vdb.athena as athena
from biograph.tools.log import debug

BIOGRAPH = './python/functest/biograph_main'
GOLDEN_DIR = 'golden/ftest'

VDB_DEFAULTS = {
    'VDB_DB': f'vdbtest{str(hex(int(time()*10000000))[8:])}',
    'VDB_BUCKET': 'vdb-test'
}

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Set VDB environment variables to the default if running under bazel not already set.
    if "TEST_WORKSPACE" in os.environ:
        for var, value in VDB_DEFAULTS.items():
            if var not in os.environ:
                print(f"Runinng under bazel; setting default value for {var}")
                os.environ[var] = value
    else:
        print("Running manually; not setting VDB environment.")
        # pylint: disable=global-statement
        global BIOGRAPH
        BIOGRAPH = 'biograph'
    print("VDB configuration:")
    for var in VDB_DEFAULTS:
        print(f"  {var}={os.environ[var]}")
    print(f"Using command '{BIOGRAPH}' to run biograph")

class VDBTestCases(unittest.TestCase):
    """ unittest Test definitions follow """

    # keep pylint and sha1_file happy
    data_dir = None

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit code is correct, return STDOUT """
        print(f"Running {cmd}", flush=True)
        proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        self.assertEqual(proc.returncode, code, f"{cmd} exited code {proc.returncode}, expected {code}")
        return proc.stdout.decode('utf-8').strip()

    # @unittest.skip('skipping vdb_basic')
    def test_000_vdb_basic(self):
        """ basic vdb commands """
        cmd = f"{BIOGRAPH} vdb"

        # exit 1
        self.check(cmd, code=1)

        # exit 0
        self.check(f"{cmd} -h")

        # other commands
        commands = {
            "vcf": ["import", "export", "list", "delete", "sort"],
            "study": [
                "create",
                # "meta",
                "add",
                "filter",
                "list",
                "show",
                "delete",
                "freeze",
                "unfreeze",
                "export",
            ]
        }
        for cat in commands:
            for sub in commands[cat]:
                self.check(f"{cmd} {cat} {sub} -h")

    # @unittest.skip('skipping vdb_study')
    def test_001_vdb_study(self):
        """ Unique study id for this test run """
        # create a study
        study_id = "vdb_study_test_" + str(hex(int(time()*10000000))[8:])

        self.check(f"{BIOGRAPH} vdb study create {study_id}")
        # verify it exists
        self.check(f"{BIOGRAPH} vdb study show {study_id}")
        # can't create it again
        self.check(f"{BIOGRAPH} vdb study create {study_id}", code=1)
        # freeze it
        self.check(f"{BIOGRAPH} vdb study freeze {study_id}")
        # can't delete it
        self.check(f"{BIOGRAPH} vdb study delete {study_id}", code=1)
        # unfreeze it
        self.check(f"{BIOGRAPH} vdb study unfreeze {study_id}")
        # delete it
        self.check(f"{BIOGRAPH} vdb study delete {study_id}")
        # it is gone
        self.check(f"{BIOGRAPH} vdb study show {study_id}", code=1)

    # @unittest.skip('skipping vdb_vcf')
    def test_002_vdb_vcf(self):
        """ vdb vcf utils """
        with tempfile.TemporaryDirectory() as tmpdir:
            # compressed vcf
            self.check(f"{BIOGRAPH} vdb vcf sort -i {GOLDEN_DIR}/vdb/vdb002.vcf.gz -o {tmpdir}/sorted.vcf")
            # uncompressed vcf on stdin and randomize the variants
            self.check(f"(zgrep ^# {GOLDEN_DIR}/vdb/vdb002.vcf.gz ; zgrep -v ^# {GOLDEN_DIR}/vdb/vdb002.vcf.gz | shuf) | {BIOGRAPH} vdb vcf sort -o {tmpdir}/sorted2.vcf")
            # should be identical
            self.assertEqual(sha1_file(f"{tmpdir}/sorted.vcf"), sha1_file(f"{tmpdir}/sorted2.vcf"))

            # uncompressed vcf, chromosomal order
            self.check(f"{BIOGRAPH} vdb vcf sort -i {tmpdir}/sorted.vcf -c -o {tmpdir}/sorted-chrom.vcf")

            # compare vs. vcf-sort
            self.check(f"zcat {GOLDEN_DIR}/vdb/vdb002.vcf.gz | vcf-sort > {tmpdir}/vcf-sorted.vcf")
            self.check(f"zcat {GOLDEN_DIR}/vdb/vdb002.vcf.gz | vcf-sort -c > {tmpdir}/vcf-sorted-chrom.vcf")

            self.assertEqual(sha1_file(f"{tmpdir}/sorted.vcf"), sha1_file(f"{tmpdir}/vcf-sorted.vcf"))
            self.assertEqual(sha1_file(f"{tmpdir}/sorted-chrom.vcf"), sha1_file(f"{tmpdir}/vcf-sorted-chrom.vcf"))

    # @unittest.skip('skipping annotations')
    def test_003_vdb_annotations(self):
        """ Annotation tests """
        cmd = f"{BIOGRAPH} vdb"
        random_id = str(hex(int(time()*10000000))[8:])

        # import some samples
        test_ids = {}

        test_ids[2] = self.check(f"{cmd} vcf import -s VDB002_{random_id} {GOLDEN_DIR}/vdb/vdb002.vcf.gz")
        test_ids[4] = self.check(f"{cmd} vcf import -s VDB004_{random_id} {GOLDEN_DIR}/vdb/vdb004.vcf.gz")

        self.check(f"{cmd} study create {random_id}")
        # wildcard add
        self.check(f"{cmd} study add {random_id} '*_{random_id}'")
        # can't add the same samples twice
        self.check(f"{cmd} study add {random_id} '*_{random_id}'", code=1)

        anno_aid = self.check(f"{cmd} anno import -i {GOLDEN_DIR}/vdb/miniclinvar.vcf.gz MiniClinVar 2021-03-25 --refhash hs37d5")

        with tempfile.TemporaryDirectory() as tmpdir:
            # no header since it contains the random study and sample IDs
            self.check(f"{cmd} study export {random_id} --no-header --anno MiniClinVar -o {tmpdir}/anno.vcf")
            # exported correctly with annotations
            self.assertEqual(sha1_file(f"{tmpdir}/anno.vcf"), "b2d639cffe095ab4c0297a2310a31ad563221d36")

        self.check(f"{cmd} anno delete {anno_aid}")
        # can only delete it once
        self.check(f"{cmd} anno delete {anno_aid}", code=1)

    # @unittest.skip('skip import_export')
    def test_004_vdb_import_export(self): # pylint: disable=too-many-statements
        """
        vdb import export

        This test also includes basic import and caching tests, to cut
        down on excessive test time.
        """
        cmd = f"{BIOGRAPH} vdb"
        random_id = str(hex(int(time()*10000000))[8:])
        study_id = f"vdb_import_export_test_{random_id}"

        test_ids = {}

        test_ids[2] = self.check(f"{cmd} vcf import -s VDB002_{random_id} {GOLDEN_DIR}/vdb/vdb002.vcf.gz")
        test_ids[3] = self.check(f"{cmd} vcf import -s VDB003_{random_id} {GOLDEN_DIR}/vdb/vdb003.vcf.gz")
        test_ids[4] = self.check(f"{cmd} vcf import -s VDB004_{random_id} {GOLDEN_DIR}/vdb/vdb004.vcf.gz")

        # this really needs to use aid instead of sample name
        self.check(f"{cmd} study create {study_id}")

        # Wildcard, but VDB003 reference doesn't match
        self.check(f"{cmd} study add {study_id} VDB*_{random_id}", code=1)

        # success
        self.check(f"{cmd} study add {study_id} VDB002_{random_id}")
        # reference doesn't match
        self.check(f"{cmd} study add {study_id} VDB003_{random_id}", code=1)
        # add by aid
        self.check(f"{cmd} study add {study_id} {test_ids[4]}")
        # can't add the same aid twice
        self.check(f"{cmd} study add {study_id} {test_ids[4]}", code=1)

        # Check VDB003 separately for correct contig style
        self.check(f"{cmd} study create {study_id}003")
        self.check(f"{cmd} study add {study_id}003 VDB003_{random_id}")

        # caching
        db = athena.connect()
        with tempfile.TemporaryDirectory() as tmpdir:
            now = time()
            query = f"""
                        SELECT sv.pos, sv.sample FROM study_variants sv
                        JOIN variants v ON
                        sv.chrom = v.chrom
                        AND sv.alt = v.alt
                        AND sv.aid = v.aid
                        WHERE sv.study_name = %(study_id)s
                    ;
                    """
            debug(query)
            # study_id is unique, so cache=True will execute and save to cache
            db.query_fetch_csv(query, f"{tmpdir}/uncached.csv", params={"study_id": study_id}, cache=True)
            delta = time() - now
            logging.warning(f"uncached fetch in {delta}")
            self.assertTrue(delta > 3, f"non-cached query finished implausibly quickly ({delta} seconds)")
            now = time()
            # fetch from cache
            db.query_fetch_csv(query, f"{tmpdir}/cached.csv", params={"study_id": study_id}, cache=True)
            delta = time() - now
            logging.warning(f"cached fetch in {delta}")
            self.assertTrue(delta < 2, f"cached query took too long ({delta} seconds)")
            # files should be identical
            self.assertEqual(sha1_file(f"{tmpdir}/uncached.csv"), sha1_file(f"{tmpdir}/cached.csv"))

        with tempfile.TemporaryDirectory() as tmpdir:
            for i in (2, 3, 4):
                self.check(f"{cmd} vcf export --aid {test_ids[i]} -o {tmpdir}/VDB00{i}.vcf.gz")
                self.assertEqual(sha1_file(f"{GOLDEN_DIR}/vdb/vdb00{i}.vcf.gz"), sha1_file(f"{tmpdir}/VDB00{i}.vcf.gz"))

            self.check(f"{cmd} study export {study_id} -o {tmpdir}/{study_id}.vcf")

            # correct distribution of N_MISS and F_MISS
            self.assertEqual(
                self.check(f"bcftools query -f '%N_MISS %F_MISS\n' {tmpdir}/{study_id}.vcf | sort | uniq -c | tr '\\n' ' '"),
                "1039 0 0    2922 1 0.5"
            )

            for f in [2, 4]:
                self.check(f"zcat {GOLDEN_DIR}/vdb/vdb00{f}.vcf.gz | grep -v '^#' | cut -f 1-5,9,10 > {tmpdir}/00{f}.orig.vcf")

            self.check(f"grep -v ^# {tmpdir}/{study_id}.vcf | cut -f 1-5,9,10 | grep -v '\\.$' > {tmpdir}/002.from_merge.vcf")
            self.check(f"grep -v ^# {tmpdir}/{study_id}.vcf | cut -f 1-5,9,11 | grep -v '\\.$' > {tmpdir}/004.from_merge.vcf")

            self.check(f"zcat {GOLDEN_DIR}/vdb/vdb003.vcf.gz | grep -v ^# | cut -f 1-5,9,10 > {tmpdir}/003.orig.vcf")
            self.check(f"{cmd} study export {study_id}003 | grep -v ^# | cut -f 1-5,9,10 > {tmpdir}/003.from_study.vcf")

            # 002 should be identical
            self.assertEqual(sha1_file(f"{tmpdir}/002.orig.vcf"), sha1_file(f"{tmpdir}/002.from_merge.vcf"))

            # 003 should be identical
            self.assertEqual(sha1_file(f"{tmpdir}/003.orig.vcf"), sha1_file(f"{tmpdir}/003.from_study.vcf"))

            # 004 is nearly identical, but includes the MP format field (and :. for sample)
            self.assertEqual(sha1_file(f"{tmpdir}/004.from_merge.vcf"), "4fa5f9c091eb64e06a6b66b2343af26facb0d228")

            # Can't square off a sample that is not in the study
            self.check(f"{cmd} study export {study_id} --no-header --square-off VDB003_{random_id}", 1)

            # 002 with and without --square-off should be identical
            self.check(f"{cmd} study export {study_id} --no-header --square-off VDB002_{random_id} -o {tmpdir}/square.vcf")
            self.check(f"grep -v ^# {tmpdir}/{study_id}.vcf | cut -f 1-6,10 > {tmpdir}/002.not.square.vcf")

            # this also implicitly tests --no-header
            self.check(f"cut -f 1-6,10 < {tmpdir}/square.vcf > {tmpdir}/002.square.vcf")
            self.assertEqual(sha1_file(f"{tmpdir}/002.square.vcf"), sha1_file(f"{tmpdir}/002.not.square.vcf"))

            # This includes a 1bp deletion at 16097167, which ends at 16097168
            self.check(f"""{cmd} study filter {study_id}003 --include "chrom = '22' && pos >= 16096307 and varend <= 16097168" """)
            self.assertEqual(self.check(f"{cmd} study export {study_id}003 --no-header | wc -l"), "3")

            # The one in the middle is lowq
            self.check(f"""{cmd} study filter {study_id}003 --exclude "FILTER != 'PASS'" """)
            self.assertEqual(self.check(f"{cmd} study export {study_id}003 --no-header | wc -l"), "2")

            # The SNP is only qual 33
            self.check(f"""{cmd} study filter {study_id}003 --include "qual > 60" """)
            self.assertEqual(self.check(f"{cmd} study export {study_id}003 --no-header | wc -l"), "1")

            # Roll back to checkpoint. No need to remerge!
            self.assertEqual(self.check(f"{cmd} study export {study_id}003 --no-header --checkpoint 3 | wc -l"), "2")

            # More range tests
            self.check(f"""{cmd} study filter {study_id} --include 'CHROM = "22" AND sample_name = "VDB002_{random_id}"' """)
            self.assertEqual(self.check(f"{cmd} study export {study_id} --no-header | wc -l"), "100")

            # Copy a study at a checkpoint
            self.check(f"""{cmd} study create {study_id}_copy""")
            self.check(f"""{cmd} study add {study_id}_copy --from {study_id} --checkpoint 2 '*'""")
            # export with --fields
            self.check(f"""{cmd} study export {study_id}_copy -o {tmpdir}/copy.vcf --fields GT:DP:AD --no-header""")

            # revert
            self.check(f"""{cmd} study revert {study_id} --checkpoint 2""")

            # --remerge required since fields changed
            self.check(f"""{cmd} study export {study_id} -o {tmpdir}/orig.vcf --fields GT:DP:AD --no-header --remerge""")
            self.assertEqual(sha1_file(f"{tmpdir}/orig.vcf"), sha1_file(f"{tmpdir}/copy.vcf"))

            # add another sample
            aid = self.check(f"""{cmd} vcf import {GOLDEN_DIR}/vdb/vdb003.hs37d5.vcf.gz""")
            # add variants from vdb003 hs37d5
            self.check(f"""{cmd} study add {study_id} {aid}""")
            # F_MISS: remove variants missing in 40% of samples
            self.check(f"""{cmd} study filter {study_id} --include 'F_MISS < 0.4'""")
            # S_MISS: eliminate all but vdb002
            self.check(f"""{cmd} study filter {study_id} --exclude 'S_MISS > 0.2'""")

            self.check(f"""{cmd} study export {study_id} -o {tmpdir}/s_miss.vcf""")
            # check variant count
            self.assertEqual(self.check(f"grep -v ^# {tmpdir}/s_miss.vcf | wc -l"), "1466")
            # should only contain VDB002
            self.assertEqual(self.check(f"grep ^#CHROM {tmpdir}/s_miss.vcf | cut -f 10-"), f"VDB002_{random_id}")

            # POS 1-to-0 translation. This matches exactly one variant on chr1 for VDB002
            self.check(f"""{cmd} study filter {study_id} --include 'POS = 10622'""")
            self.assertEqual(
                self.check(f"""{cmd} study export {study_id} --no-header | cut -f 1-5"""),
                "1\t10622\t.\tTT\tGCGCAGGC"
            )

    # @unittest.skip('skip cleanup')
    def test_999_cleanup(self):
        ''' delete the database and s3 data '''
        if not os.environ["VDB_DB"].startswith("vdbtest"):
            self.assertEqual(f"Cowardly refusing to delete test database {os.environ['VDB_DB']}", "")

        db = athena.connect()
        futures = []
        for table_class in [vars(tc).values() for tc in vars(db.table).values()]:
            for table_name in table_class:
                print(f"DROP TABLE {table_name};")
                futures.append(db.query(f"DROP TABLE {table_name};", block=False))

        for future in futures:
            future.result()

        db.query(f"DROP DATABASE {db.database} CASCADE;")

        # s3 is last so results/ gets removed after dropping the database
        db.s3_rm_recursive(db.database)

if __name__ == '__main__':
    unittest.main(verbosity=2)
