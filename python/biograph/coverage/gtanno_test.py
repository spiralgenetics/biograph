# pylint: disable=missing-docstring
import unittest
import biograph
from biograph.coverage import GTAnno, VcfEntryInfo
import biograph.variants as bgexvar

class CovAnnoTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/proband_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")
        cls.vcf_file = "datasets/lambdaToyData/benchmark/family.vcf.gz"
        cls.asms = VcfEntryInfo.parse_region(cls.vcf_file, 'lambda', 0,
                                             cls.ref.scaffold_lens["lambda"],
                                             sample_index=None)
        # terribly verbose answer checking
        cls.answer_vcf = []
        cls.answer_fmt = []
        cls.answer_vcf.append(['lambda', '25054', '.', 'G', 'T', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:28475:144:147:0,147:107:0,107', '.', '.'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 5, 'PL': '99,72,5', 'DP': '120', 'AD': '0,120', 'RC': '0'})
        cls.answer_vcf.append(['lambda', '30278', '.', 'G',
                               'GGCCTAGGCGGGAACGTGGGCCATGGTGGCTGCCGCATGTACTGGCGATTGATCCTCCTGCAACCTGAAGGGACGGCCGCGGGAACGTCTCCGATAATGAAGGCTTGCACTCATATACTATCCAAGCCACGGGTGATACACCCGTGGCACTAAGAACGTTATAGAGAACCTATCTTTCGGGGATGGGCCTATTGCGTCTAACATAGACACTTTAAGGCTAATGAAGTTTGTAGCTAAGACCGCTGGGGAGTGAATAGCGGGACACGAATGGTCGGGAAGCAAAACGAAACGGAGGATTCTC',
                               '100', 'PASS', 'NS=1;END=30278;SVLEN=300;SVTYPE=INS', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:36168:142:119:0,119:35:0,35', '.', '1/1:1|1:100:36048:145:124:0,124:42:0,42'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 10, 'PL': '99,99,10', 'DP': '224', 'AD': '0,224', 'RC': '2'})
        cls.answer_vcf.append(['lambda', '33538', '.', 'A', 'AAG', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:40313:143:135:0,135:100:0,100', '1/1:1|1:100:41771:146:158:0,158:122:0,122', '.'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 4, 'PL': '99,66,4', 'DP': '117', 'AD': '1,116', 'RC':'0'})
        cls.answer_vcf.append(['lambda', '37744', '.', 'AG', 'A', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:12941:145:146:0,146:102:0,102', '.', '1/1:1|1:100:45276:144:163:0,163:125:0,125'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 5, 'PL': '99,67,5', 'DP': '112', 'AD': '0,112', 'RC': '0'})
        cls.answer_vcf.append(['lambda', '41165', '.', 'GAC', 'G', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:8704:146:168:0,168:123:0,123', '1/1:1|1:100:51266:144:146:0,146:116:0,116', '.'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 6, 'PL': '99,78,6', 'DP': '129', 'AD': '0,129', 'RC': '0'})
        cls.answer_vcf.append(['lambda', '43279', '.', 'G', 'C', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:6241:145:117:0,117:88:0,88', '1/1:1|1:100:53789:145:146:0,146:113:0,113', '.'])
        cls.answer_fmt.append({'GT': '1/1', 'GQ': 4, 'PL': '99,55,4', 'DP': '91', 'AD': '0,91', 'RC':'0'})

    def test_lambda(self):
        anno = GTAnno()
        asms = bgexvar.add_ref_assemblies(self.ref, 'lambda', self.asms, 1300)

        asms = bgexvar.generate_read_cov(self.rm, asms)
        asms = bgexvar.generate_pair_cov(self.rm, asms)
        data = anno.parse(self.rm, asms)
        found = 0
        trimmed_answer = [fields[:9] for fields in self.answer_vcf]
        for i in data:
            if i.vcf_entry_info.vcf_entry in trimmed_answer:
                idx = trimmed_answer.index(i.vcf_entry_info.vcf_entry)
                self.assertEqual(self.answer_fmt[idx], i.vcf_entry_info.new_fmt)
                found += 1
        self.assertEqual(found, len(self.answer_vcf))


if __name__ == '__main__':
    unittest.main(verbosity=2)
