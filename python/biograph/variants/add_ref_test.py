# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

class ReadCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")
        cls.reflen = cls.ref.scaffold_lens['lambda']

    def test_add_ref(self):
        asms = bgexvar.add_ref_assemblies(self.ref, "lambda", [], whole_ref=True, max_len=100)
        asms = list(asms)
        self.assertEqual(len(asms), (self.reflen // 100) + 1)

        seq = biograph.Sequence()
        for a in asms:
            seq += a.seq
        self.assertEqual(len(seq), self.reflen)

        rc_asms = bgexvar.add_ref_assemblies(self.ref, "lambda", [], whole_ref=True, max_len=100,
                                             rev_comp=True)
        rc_asms = list(rc_asms)
        self.assertEqual(len(rc_asms), (self.reflen // 100) + 1)
        rc_seq = biograph.Sequence()
        for a in rc_asms:
            rc_seq += a.seq
        self.assertEqual(len(rc_seq), self.reflen)

        self.assertEqual(seq, rc_seq.rev_comp())

if __name__ == '__main__':
    unittest.main(verbosity=2)
