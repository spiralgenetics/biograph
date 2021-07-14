"""
basic.py: Pre-run sanity checks.
"""

from __future__ import print_function

import unittest
import subprocess

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
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

    def check(self, cmd, code=0):
        """ Run a shell command, assert exit is zero. """
        self.assertEqual(code, subprocess.call(cmd + " >/dev/null 2>&1", shell=True), "%s exited non-zero" % cmd)

    def test_biograph(self):
        """
            biograph
        """
        self.check('bgbinary', 1)
        self.check('bgbinary blastn', 1)

if __name__ == '__main__':
    unittest.main(verbosity=2)
