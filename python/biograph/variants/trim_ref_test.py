# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

class AssemblyTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def test_insert(self):
        asms = [bgexvar.Assembly(2666, 2667, "CA", 123)]
        trimmed = list(bgexvar.trim_ref(self.ref, "lambda", asms))
        self.assertEqual(len(trimmed), 1)
        a = trimmed[0]
        self.assertEqual(a.assembly_id, 123)
        self.assertEqual(a.left_offset, 2667)
        self.assertEqual(a.right_offset, 2667)
        self.assertEqual(a.seq, "A")

if __name__ == '__main__':
    unittest.main(verbosity=2)
