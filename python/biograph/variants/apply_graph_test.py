# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

class ApplyGraphTracker: # pylint: disable=too-few-public-methods
    def __init__(self):
        self.invocations = []

    def on_graph_context(self, ctx):
        assert(not ctx.a.matches_reference)
        inv = (
            str(ctx.left_ref.seq),
            str(ctx.a.seq),
            str(ctx.right_ref.seq),
            [str(ref.seq) for ref in ctx.refs],
            len(ctx.edge_coverage(ctx.a.read_coverage, ctx.ref_coverage()).variant_start),
        )
        self.invocations.append(inv)


def vcf_assembly(pos, ref, alt, asm_id):
    pos = int(pos)-1
    if ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)

class ApplyGraphTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def test_apply_graph(self):
        self.maxDiff = None
        asms = [vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1),
                vcf_assembly("2667", "C", "CA", 2)]
        pc = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 105)
        pc = list(bgexvar.generate_read_cov(self.rm, pc))
        orig_pc_len = len(pc)
        seen = ApplyGraphTracker()
        pc = list(bgexvar.apply_graph(pc, seen.on_graph_context))
        self.assertEqual(len(pc), orig_pc_len)

        self.assertEqual([
            # delete at 2191
            ('ATCCCGTATCTGCTGGGATACTGGCGGGATTGACCCGACCATTGTGTATGAACGCTCGAAAAAACATGGGCTGTTCCGGGTGATCCCCATTAAAGGGGCATCCGT',
             '',
             'CGCTTCACACTGACGCCGGAAGGGGATGAACCGCTTCCCGGTGCCGTTCACTTCCCGAATAACCCGGATATTTTTGATCTGACCGAAGCGCAGCAGCTGACTGCT',
             ['CTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC'],
             70,
            ),
            # insert at 2667
            ('AAGAGGATGGTGCAGCAACCAACAAGAAAACACTGGCAGATTACGCCCGTGCCTTATCCGGAGAGGATGAATGACGCGACAGGAAGAACTTGCCGCTGCCCGTGC',
             'A',
             'GGCACTGCATGACCTGATGACAGGTAAACGGGTGGCAACAGTACAGAAAGACGGACGAAGGGTGGAGTTTACGGCCACTTCCGTGTCTGACCTGAAAAAATATAT',
             [],
             69,
            ),
        ], seen.invocations)

if __name__ == '__main__':
    unittest.main(verbosity=2)
