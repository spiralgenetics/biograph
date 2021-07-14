"""
sam.py: SAM file parser
"""
from __future__ import print_function

import unittest

from python.functest.utils.sam import (
    sam_format,
    parse_sam_line
)

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)

from python.functest.utils.defaults import GOLDEN_DIR

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class UtilitiesTestCases(unittest.TestCase):
    """    unittest Test definitions follow """

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def test_sam_parser_success(self):
        """
            human_100.sam should have 100 entries in it, and parsing should throw no exceptions.
        """
        # self.cleanup = False

        sam_file = '%s/../human_100.sam' % GOLDEN_DIR
        count = 0

        with open(sam_file, 'r') as sam:
            for line in sam:
                # Skip tags
                if line[0] == '@':
                    continue

                alignment = parse_sam_line(line)

                # Verify that the type conversions are all correct
                types = {}
                for entry in sam_format():
                    types[entry['name']] = entry['type']

                for field in alignment:
                    self.assertIs(type(alignment[field]), types[field])

                count = count + 1

        self.assertEqual(count, 100)


    def test_sam_parser_throw(self):
        """
            Parsing a non-SAM file should throw an exception.
        """
        # self.cleanup = False

        some_file = '%s/fake_results' % GOLDEN_DIR

        try:
            with open(some_file, 'r') as something:
                for line in something:
                    parse_sam_line(line)
        # pylint: disable=broad-except
        except Exception:
            pass
        else:
            self.fail('Exception should have been called when parsing a non-SAM file.')


if __name__ == '__main__':
    unittest.main(verbosity=2)
