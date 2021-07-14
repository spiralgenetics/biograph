# pylint: disable=missing-docstring
import unittest
import biograph
from biograph.coverage import CovAnno, VcfEntryInfo
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
                                             cls.ref.scaffold_lens["lambda"], sample_index=None)
        # terribly verbose answer checking
        cls.answer_vcf = []
        cls.answer_fmt = []
        cls.answer_vcf.append(['lambda', '2191', '.',
                               'TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAAC', 'T',
                               '100', 'PASS', 'NS=1;END=2291;SVLEN=100;SVTYPE=DEL', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:54842:144:126:0,126:97:0,97', '.', '1/1:1|1:100:57915:145:143:0,143:115:0,115'])
        cls.answer_fmt.append({'US': '2,74', 'DS': '0,74', 'UC': '2,127', 'DC': '0,126',
                               'UDC': '1,57', 'UCC': '1,70', 'DDC': '0,56', 'DCC': '0,70',
                               'UMO': '149,141', 'DMO': '0,141', 'UXO': '149,150', 'DXO':
                               '0,150', 'NR': '2,127', 'MO': '0,141', 'XO': '0,150',
                               'XC': '0,127', 'AC': '0,126', 'EC': '0,126', 'MC': '0,126', 'MP': '2,2'})
        cls.answer_vcf.append(['lambda', '2667', '.', 'C', 'CA', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:54495:146:153:0,153:115:0,115', '.', '1/1:1|1:100:2062:145:152:0,152:119:0,119'])
        cls.answer_fmt.append({'US': '1,75', 'DS': '0,74', 'UC': '2,155', 'DC': '0,155',
                               'UDC': '0,88', 'UCC': '2,67', 'DDC': '0,88', 'DCC': '0,67',
                               'UMO': '150,146', 'DMO': '0,146', 'UXO': '150,150', 'DXO':
                               '0,150', 'NR': '2,157', 'MO': '0,146', 'XO': '0,150',
                               'XC': '0,155', 'AC': '0,154', 'EC': '0,155', 'MC': '0,153', 'MP': '2,2'})
        cls.answer_vcf.append(['lambda', '5897', '.', 'G', 'A', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:6698:145:160:0,160:118:0,118', '.', '1/1:1|1:100:53723:146:153:0,153:105:0,105'])
        cls.answer_fmt.append({'US': '0,74', 'DS': '0,75', 'UC': '0,160', 'DC': '0,161',
                               'UDC': '0,76', 'UCC': '0,84', 'DDC': '0,76', 'DCC': '0,85',
                               'UMO': '0,145', 'DMO': '0,145', 'UXO': '0,150', 'DXO':
                               '0,150', 'NR': '0,161', 'MO': '0,145', 'XO': '0,150', 'XC':
                               '0,161', 'AC': '0,160', 'EC': '0,160', 'MC': '0,160', 'MP': '2,2'})
        cls.answer_vcf.append(['lambda', '7146', '.', 'G', 'GTA', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:49155:144:144:0,144:110:0,110', '.', '1/1:1|1:100:52298:144:138:0,138:99:0,99'])
        cls.answer_fmt.append({'US': '2,75', 'DS': '1,72', 'UC': '4,144', 'DC': '2,148',
                               'UDC': '1,71', 'UCC': '3,73', 'DDC': '0,73', 'DCC': '2,75',
                               'UMO': '149,140', 'DMO': '150,140', 'UXO': '150,150', 'DXO':
                               '150,150', 'NR': '6,152', 'MO': '150,140', 'XO': '150,150',
                               'XC': '2,148', 'AC': '2,146', 'EC': '2,145', 'MC': '2,144', 'MP': '2,2'})
        cls.answer_vcf.append(['lambda', '8493', '.', 'C', 'CT', '100', 'PASS', 'NS=1', 'GT:PG:GQ:PI:OV:DP:AD:PDP:PAD',
                               '1/1:1|1:100:47631:146:171:0,171:135:0,135', '1/1:1|1:100:48926:145:139:0,139:106:0,106', '.'])
        cls.answer_fmt.append({'US': '1,75', 'DS': '0,73', 'UC': '2,172', 'DC': '0,171',
                               'UDC': '1,88', 'UCC': '1,84', 'DDC': '0,89', 'DCC': '0,82',
                               'UMO': '150,145', 'DMO': '0,145', 'UXO': '150,150', 'DXO':
                               '0,150', 'NR': '2,174', 'MO': '0,145', 'XO': '0,150',
                               'XC': '0,172', 'AC': '0,171', 'EC': '0,171', 'MC': '0,171', 'MP': '2,2'})

    def test_lambda(self):
        self.maxDiff = 10000
        anno = CovAnno(True)
        asms = bgexvar.add_ref_assemblies(self.ref, 'lambda', self.asms, 250)
        data = anno.parse(self.rm, asms)
        found = 0
        expecteds, actuals = [], []
        trimmed_answer = [fields[:9] for fields in self.answer_vcf]
        for i in data:
            if i.matches_reference:
                continue
            if i.vcf_entry_info.vcf_entry in trimmed_answer:
                idx = trimmed_answer.index(i.vcf_entry_info.vcf_entry)

                expecteds.append(self.answer_fmt[idx])
                actuals.append(i.vcf_entry_info.new_fmt)
                print(self.answer_fmt[idx])
                found += 1
        self.assertEqual(expecteds, actuals)

        self.assertEqual(found, len(self.answer_vcf))

if __name__ == '__main__':
    unittest.main(verbosity=2)
