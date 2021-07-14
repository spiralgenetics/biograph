# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import re

import biograph
import biograph.variants as bgexvar
from biograph.coverage import SamOutput

def vcf_assembly(pos, ref, alt, asm_id):
    pos = int(pos)-1
    if ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)

CIGAR_RE = re.compile(r'(\d+)(\D+)')
def parse_cigar(cigarstring):
    cigars = CIGAR_RE.split(cigarstring)
    if (not cigars) or cigars[0]:
        raise RuntimeError(f"Could not parse cigar {cigarstring}")
    cigars = cigars[1:]
    while cigars:
        bases, op, between = cigars[:3]
        if between:
            raise RuntimeError(f"Could not parse cigar {cigarstring}")
        yield op, int(bases)
        cigars = cigars[3:]

class AlignReadsCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/mother_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")
        cls.short_asm_list = [
            vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1001),
            vcf_assembly("2667", "C", "CA", 1002),
            vcf_assembly("20056", "C", "CAAGGGCAACTCGCCGTACATCTATGGGACATTTTCCGCGACCGTTCGGAACTTTCGTCTGCCGAGTCAGGATAAAGGCTGATTGAGTTTGGGAGAAAGCA", 1009)
        ]
        cls.long_asm_list = [
            vcf_assembly("2191", "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC", "T", 1001),
            vcf_assembly("2667", "C", "CA", 1002),
            vcf_assembly("5897", "G", "A", 1003),
            vcf_assembly("7146", "G", "GTA", 1004),
            vcf_assembly("8493", "C", "CT", 1005),
            vcf_assembly("9495", "C", "G", 1006),
            vcf_assembly("13383", "T", "A", 1007),
            vcf_assembly("16619", "ACCGCAGTAATGTGGTGATGCCGGATGATGGCGCGCCGTTCCGCTACAGCTTCAGCGCCCTGAAGGACCGCCATAATGCCGTTGAGGTGAACTGGATTGACCCGAACAACGGCTGGGAGACGGCGACAGAGCTTGTTGAAGATACGCAGGCCATTGCCCGTTACGGTCGTAATGTTACGAAGATGGATGCCTTTGGCTGTACCAGCCGGGGGCAGGCACACCGCGCCGGGCTGTGGCTGATTAAAACAGAACTGCTGGAAACGCAGACCGTGGATTTCAGCGTCGGCGCAGAAGGGCTTCG", "A", 1008),
            vcf_assembly("20056", "C", "CAAGGGCAACTCGCCGTACATCTATGGGACATTTTCCGCGACCGTTCGGAACTTTCGTCTGCCGAGTCAGGATAAAGGCTGATTGAGTTTGGGAGAAAGCA", 1009),
            vcf_assembly("24304", "TGA", "T", 1010),
            vcf_assembly("25054", "G", "T", 1011),
            vcf_assembly("30278", "G", "GGCCTAGGCGGGAACGTGGGCCATGGTGGCTGCCGCATGTACTGGCGATTGATCCTCCTGCAACCTGAAGGGACGGCCGCGGGAACGTCTCCGATAATGAAGGCTTGCACTCATATACTATCCAAGCCACGGGTGATACACCCGTGGCACTAAGAACGTTATAGAGAACCTATCTTTCGGGGATGGGCCTATTGCGTCTAACATAGACACTTTAAGGCTAATGAAGTTTGTAGCTAAGACCGCTGGGGAGTGAATAGCGGGACACGAATGGTCGGGAAGCAAAACGAAACGGAGGATTCTC", 1012),
            vcf_assembly("33538", "A", "AAG", 1013),
            vcf_assembly("37744", "AG", "A", 1014),
            vcf_assembly("41165", "GAC", "G", 1015),
            vcf_assembly("43279", "G", "C", 1016),
            vcf_assembly("45198", "GA", "G", 1017)
        ]

    @staticmethod
    def consume_asms(asms):
        for _ in asms:
            pass

    def test_sequences(self):
        asms = self.long_asm_list
        sam = SamOutput(self.rm)
        sam.add_reference(self.ref)

        asms = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 200)
        asms = bgexvar.generate_read_cov(self.rm, asms)

        def on_aligned_single(read_id, aligned):
            read_seq = self.rm.get_read_by_id(read_id).get_seqset_entry().sequence()
            self.assertEqual(aligned.seq, read_seq)
            seq_len = 0
            ref_len = 0
            cigar = parse_cigar(aligned.cigar)
            for cigar_op, bases in cigar:
                self.assertGreater(bases, 0)
                if cigar_op == "M":
                    seq_len += bases
                    ref_len += bases
                elif cigar_op in ('D', 'N'):
                    ref_len += bases
                elif cigar_op == "I":
                    seq_len += bases
                else:
                    self.assertEqual(cigar_op, "P")
            self.assertEqual((seq_len, ref_len),
                             (len(aligned.seq), (aligned.right_offset - aligned.left_offset)))

        def on_aligned(read_ids, aligned):
            for read_id in read_ids:
                on_aligned_single(read_id, aligned)

        asms = bgexvar.align_reads(asms, on_aligned=on_aligned)
        self.consume_asms(asms)


if __name__ == '__main__':
    unittest.main(verbosity=2)
