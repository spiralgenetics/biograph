# pylint: disable=missing-docstring

from __future__ import print_function

from operator import attrgetter
import unittest
import biograph
import biograph.variants as bgexvar

class ApplyEdgesTracker:
    def __init__(self):
        self.invocations = []

    def on_assembly_edges(self, offset, left_edges, inserts, right_edges):
        self.invocations.append((offset,
                                 self.track_asms(left_edges),
                                 self.track_asms(inserts),
                                 self.track_asms(right_edges)))

    @staticmethod
    def track_asms(asms):
        trackinfo = []
        asms.sort(key=attrgetter("left_offset"))
        asms.sort(key=attrgetter("right_offset"))
        asms.sort(key=attrgetter("matches_reference"))
        for asm in asms:
            if asm.matches_reference:
                asm_ti = ("REF", asm.left_offset, asm.right_offset)
            else:
                asm_ti = ("VAR", asm.left_offset, asm.right_offset)
            trackinfo.append(asm_ti)
        return trackinfo

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

    def test_apply_edges(self):
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1),
                vcf_assembly("2667", "C", "CA", 2)]
        pc = list(bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 10))
        orig_pc_len = len(pc)
        edges_seen = ApplyEdgesTracker()
        pc = list(bgexvar.apply_edges(pc, edges_seen.on_assembly_edges))
        self.assertEqual(len(pc), orig_pc_len)

        self.assertEqual(edges_seen.invocations, [
            (2181,
             [],
             [],
             [('REF', 2181, 2191)]),
            (2191,
             [('REF', 2181, 2191)],
             [],
             [('VAR', 2191, 2291),
              ('REF', 2191, 2291)]),
            (2291,
             [('VAR', 2191, 2291),
              ('REF', 2191, 2291)],
             [],
             [('REF', 2291, 2301)]),
            (2301,
             [('REF', 2291, 2301)],
             [],
             []),
            (2657,
             [],
             [],
             [('REF',
               2657,
               2667)]),
            (2667,
             [('REF', 2657, 2667)],
             [('VAR', 2667, 2667)],
             [('REF', 2667, 2677)]),
            (2677,
             [('REF', 2667, 2677)],
             [],
             [])
        ])

    def test_exceptions(self):
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1),
                vcf_assembly("2667", "C", "CA", 2)]
        pc = list(bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 10))
        did_invoke = False
        def on_invoke(offset, left, inserts, right):
            nonlocal did_invoke
            did_invoke = True
            raise RuntimeError("Should be propagated")
        with self.assertRaises(RuntimeError):
            pc = list(bgexvar.apply_edges(pc, on_invoke))
            edges_seen = ApplyEdgesTracker()
            pc = list(bgexvar.apply_edges(pc, edges_seen.on_assembly_edges))
        self.assertTrue(did_invoke)

if __name__ == '__main__':
    unittest.main(verbosity=2)
