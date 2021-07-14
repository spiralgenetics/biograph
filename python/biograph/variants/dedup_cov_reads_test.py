# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

class PairEdgeCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def run_pair_edge_cov(self, asms, min_dist=0, max_dist=10000):
        pc = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, max_dist)
        pc = bgexvar.generate_read_cov(self.rm, pc)
        pc = bgexvar.generate_pair_cov(self.rm, pc, min_dist, max_dist)
        pc = bgexvar.generate_pair_edge_cov(pc)
        return list(pc)

    def test_multiple(self):
        asms = [bgexvar.Assembly(20056, 20056, "AAGGGCAACTCGCCGTACATCTATGGGACATTTTCCGCGACCGTTCGGAACTTTCGTCTGCCGAGTCAGGATAAAGGCTGATTGAGTTTGGGAGAAAGCA", 1)]
        pc = self.run_pair_edge_cov(asms)
        print(pc)
        pc_filtered = [a for a in pc if not a.matches_reference]
        self.assertEqual(len(pc_filtered), 1)
        a = pc_filtered[0]

        self.assertTrue(6481 in a.edge_coverage.variant_start.expand_to_list())

        dd = list(bgexvar.dedup_cov_reads(self.ref, "lambda", pc))
        dd_filtered = [a for a in dd if not a.matches_reference]
        self.assertEqual(len(dd_filtered), 1)

        a = dd_filtered[0]
        self.assertFalse(6481 in a.edge_coverage.variant_start.expand_to_list())

if __name__ == '__main__':
    unittest.main(verbosity=2)
