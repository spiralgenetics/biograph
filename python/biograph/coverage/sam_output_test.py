# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import tempfile
import multiprocessing

import biograph
import biograph.variants as bgexvar
from biograph.coverage import SamOutput
import pysam

def vcf_assembly(pos, ref, alt, asm_id):
    pos = int(pos)-1
    if ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)

class SamOutputCases(unittest.TestCase):
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


    def test_sam_output(self):
        asms = self.short_asm_list
        with tempfile.TemporaryDirectory() as tmpdir:
            sam_filename = tmpdir + "/output.sam"

            outqueue = multiprocessing.Queue()
            sam = SamOutput(self.rm)
            sam.add_reference(self.ref)

            with open(sam_filename, "wb") as samf:
                for line in sam.get_header():
                    samf.write(line)
                asms = bgexvar.add_ref_assemblies(self.ref, "lambda", asms, 140)
                asms = bgexvar.generate_read_cov(self.rm, asms)
                asms = sam.output_aligned(outqueue, "lambda", asms)
                for _ in asms:
                    pass
                outqueue.put(None)
                while True:
                    line = outqueue.get()
                    if line is None:
                        break
                    samf.write(line)

            actuals = []
            # Sampling of some of the 635 expecteds:
            expecteds = [
                [0, 2052, 145, [(0, 150)], 'CCTATACCCGCCGGAATGGTGCAGAAATGTCGATATCCCGTATCTGCTGGGATACTGGCGGGATTGACCCGACCATTGTGTATGAACGCTCGAAAAAACATGGGCTGTTCCGGGTGATCCCCATTAAAGGGGCATCCGTCTACGGAAAGC'],
                [0, 2134, 65, [(0, 116)], 'TGAACGCTCGAAAAAACATGGGCTGTTCCGGGTGATCCCCATTAAAGGGGCATCCGTCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTA'],
                [0, 19918, 81, [(0, 138), (1, 12), (6, 88)], 'CTGTGCCATGACGGAGGATGATGCCCGGCCGGAGGTGCTGCGTCGTCTTGAACTGATGGTGGAAGAGGTGGCGCGTAACGCGTCCGTGGTGGCACAGAGTACGGCAGACGCGAAGAAATCAGCCGGCGATGCCAGTGCAAGGGCAACTCG'],
                [0, 20025, 65, [(0, 31), (1, 100), (0, 19)], 'ACGCGAAGAAATCAGCCGGCGATGCCAGTGCAAGGGCAACTCGCCGTACATCTATGGGACATTTTCCGCGACCGTTCGGAACTTTCGTCTGCCGAGTCAGGATAAAGGCTGATTGAGTTTGGGAGAAAGCAATCAGCTGCTCAGGTCGCG'],
                [0, 20056, 145, [(6, 83), (1, 17), (0, 133)], 'GAGTTTGGGAGAAAGCAATCAGCTGCTCAGGTCGCGGCCCTTGTGACTGATGCAACTGACTCAGCACGCGCCGCCAGCACGTCCGCCGGACAGGCTGCATCGTCAGCTCAGGAAGCGTCCTCCGGCGCAGAAGCGGCATCAGCAAAGGCC']
            ]

            with pysam.AlignmentFile(sam_filename, "rb") as samf:
                for a in samf:
                    actual = [a.reference_id, a.reference_start, a.flag,
                              a.cigar, a.query_sequence]
                    actuals.append(actual)

            for expected in expecteds:
                self.assertIn(expected, actuals)
            self.assertEqual(len(actuals), 635)


if __name__ == '__main__':
    unittest.main(verbosity=2)
