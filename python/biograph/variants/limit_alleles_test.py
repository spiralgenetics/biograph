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
    a = bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)
    a.limited = False
    return a

class PhaseTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    @staticmethod
    def sort_by_asm_id(asms):
        print(f"Sorting {asms}")
        return sorted(asms, key=lambda a: a.assembly_id)

    @staticmethod
    def mark_limited(a):
        print(f"Marking limited: {a}")
        a.limited = True

    def limited_by_id(self, asms):
        limiteds = []
        for a in asms:
            while len(limiteds) <= a.assembly_id:
                limiteds.append(None)
            self.assertTrue(limiteds[a.assembly_id] is None)
            limiteds[a.assembly_id] = a.limited
        return limiteds

    def test_lambda(self):
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 0),
                vcf_assembly("2667", "C", "G", 1),
                vcf_assembly("2667", "C", "T", 2)]

        asms = bgexvar.limit_alleles(1, self.sort_by_asm_id, self.mark_limited, asms)
        asms = list(asms)
        self.assertEqual(self.limited_by_id(asms), [False, False, True])

if __name__ == '__main__':
    unittest.main(verbosity=2)
