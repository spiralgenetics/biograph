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

class ReadCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def run_read_cov(self, asms):
        pc = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 300)
        pc = bgexvar.generate_read_cov(self.rm, pc)
        return list(pc)

    def test_lambda(self):
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1),
                vcf_assembly("2667", "C", "CA", 2)]
        pc = self.run_read_cov(asms)
        self.assertEqual(len(pc), 6)
        for idx in [0, 2, 3, 5]:
            self.assertTrue(pc[idx].matches_reference)
            self.assertEqual(pc[idx].generated_by, "ADD_REF")
        self.assertEqual(pc[1].assembly_id, 1)
        self.assertEqual(pc[4].assembly_id, 2)

        self.assertEqual(pc[1].read_coverage.get_tot_read_count(), 143)
        self.assertFalse(bgexvar.ReadCoverageRead(-13, 5109, 1234) in pc[1].read_coverage)
        self.assertTrue(bgexvar.ReadCoverageRead(-13, 5109, 150) in pc[1].read_coverage)

        self.assertEqual(pc[4].read_coverage.assembly_len(), 1)
        dp = pc[4].read_coverage.calc_depths()
        self.assertEqual(dp, [152, 152])

        dp = pc[4].read_coverage.calc_depths(interbase=True)
        self.assertEqual(dp, [152, 152])
        dp = pc[4].read_coverage.calc_depths(interbase=True, include_fwd=False, readmap=self.rm)
        self.assertEqual(dp, [77, 77])
        dp = pc[4].read_coverage.calc_depths(interbase=True, include_rev=False, readmap=self.rm)
        self.assertEqual(dp, [75, 75])

        self.assertEqual(pc[4].read_coverage.get_tot_read_count(), 153)
        dp = pc[4].read_coverage.calc_depths(interbase=False)
        self.assertEqual(dp, [153])
        dp = pc[4].read_coverage.calc_depths(interbase=False, include_fwd=False, readmap=self.rm)
        self.assertEqual(dp, [77])
        dp = pc[4].read_coverage.calc_depths(interbase=False, include_rev=False, readmap=self.rm)
        self.assertEqual(dp, [76])

        with self.assertRaises(RuntimeError):
            pc[4].read_coverage.calc_depths(interbase=False, include_rev=False)

if __name__ == '__main__':
    unittest.main(verbosity=2)
