#!/usr/bin/env python3
"""
biograph_wrapper_test.py: Test the biograph wrapper commands.

This test does not (yet) run under Bazel, so run it directly from the spiral/ root
from inside a venv with biograph installed.
"""
import unittest
import subprocess
import tempfile
import os
import biograph

from pathlib import Path

GOLDEN_DIR = './golden/ftest'

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

class WrapperTestCases(unittest.TestCase):
    """    unittest Test definitions follow """

    # keep pylint happy
    data_dir = None

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit code is correct """
        actual = subprocess.call(cmd, shell=True)
        self.assertEqual(code, actual, "{0} exited code {1}".format(cmd, actual))

    # @unittest.skip('no basic')
    def test_000_biograph_basic(self):
        """ basic biograph commands """
        bgcmd = "biograph"

        # exit 1
        self.check(bgcmd, code=1)

        # exit 0
        self.check(f"{bgcmd} -h")

        # other commands
        for cmd in (
            "full_pipeline",
            "reference",
            "create",
            "discovery",
            "coverage",
            "qual_classifier",
            "vdb",
            "stats",
            "version",
            "refhash"
        ):
            self.check(f"{bgcmd} {cmd} -h")

    # @unittest.skip('nope')
    def test_001_biograph_full_pipeline(self):
        """
            biograph full_pipeline, beginning to end
        """
        bgcmd = "biograph"
        version = biograph.version().rpartition('-dev')[0]
        model = f"/share/releases/biograph_model-{version}.ml"

        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                bgcmd,
                "full_pipeline",
                f"--biograph {tmpdir}/out.bg",
                "--ref /reference/e_coli_k12_ASM584v1/",
                f"--model {model}",
                "--reads /share/bams/e_coli/e_coli_test.bam",
                "--force"
            )
            self.check(' '.join(cmd), 0)

    # @unittest.skip('nope')
    def test_002_biograph_full_pipeline_stdin(self):
        """
            biograph full_pipeline with reads on STDIN
            No need for a model if we stop at discovery
        """
        bgcmd = "biograph"

        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                "cat", "/share/bams/e_coli/e_coli_test.bam", "|",
                bgcmd,
                "full_pipeline",
                f"--biograph {tmpdir}/out.bg",
                "--ref /reference/e_coli_k12_ASM584v1/",
                "--reads -",
                "--stop", "discovery",
                "--force"
            )
            self.check(' '.join(cmd), 0)

if __name__ == '__main__':
    unittest.main(verbosity=2)
