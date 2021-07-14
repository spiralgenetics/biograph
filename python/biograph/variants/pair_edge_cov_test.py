# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar


def vcf_assembly(pos, ref, alt, asm_id):
    pos = int(pos)-1
    if ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)

class PairEdgeCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def run_pair_edge_cov(self, asms, min_dist=0, max_dist=10000):
        pc = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, max_dist*2)
        pc = bgexvar.generate_read_cov(self.rm, pc)
        pc = bgexvar.generate_pair_cov(self.rm, pc, min_dist, max_dist)
        pc = bgexvar.generate_pair_edge_cov(pc)
        pc = bgexvar.trim_ref(self.ref, "lambda", pc)
        return list(pc)

    def test_deletion(self):
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1),
                vcf_assembly("2667", "C", "CA", 2)]
        pc = self.run_pair_edge_cov(asms)
        self.assertEqual(len(pc), 2)
        self.assertEqual(pc[0].assembly_id, 1)
        self.assertEqual(pc[1].assembly_id, 2)

        self.assertEqual(len(pc[0].edge_coverage.variant_start.expand_to_list()), 121)
        self.assertTrue(5109 in list(pc[0].edge_coverage.variant_start.expand_to_list()))

    def test_obj_preserved(self):
        asms = [vcf_assembly("2667", "C", "CA", 2)]
        asms[0].random_property = "random_value"
        pc = self.run_pair_edge_cov(asms)
        self.assertEqual(len(pc), 1)
        self.assertEqual(pc[0].random_property, "random_value")
        self.assertEqual(pc[0], asms[0])

    def test_assembly_order(self):
        # Make sure out of order assemblies raise exception
        asms = [vcf_assembly("2667", "C", "CA", 2),
                vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1)]
        with self.assertRaisesRegex(RuntimeError, "Assemblies must be sorted in order"):
            self.run_pair_edge_cov(asms)

    def test_bad_assembly(self):
        # Make sure bad assemblies raise exceptions
        asm = vcf_assembly("2667", "C", "CA", 2)
        asm.left_offset = asm.right_offset + 1
        with self.assertRaisesRegex(RuntimeError, "Right offset must occur after the left offset"):
            self.run_pair_edge_cov([asm])

if __name__ == '__main__':
    unittest.main(verbosity=2)
