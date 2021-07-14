"""
filter_test.py: vdb filter parsing tests
"""
import unittest
from pyparsing import ParseException

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)

from biograph.vdb.filter import parser

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class FilterTestCases(unittest.TestCase):
    """ biograph filter tests follow """

    # keep pylint happy
    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def check_pattern(self, ts, equiv, parser_type='vdb'):
        ''' Attempt to parse the test string, then assert that the parsed output is equal to equiv '''
        result = parser(ts, parser_type)
        self.assertEqual(str(result), str(equiv))

    def test_patterns(self): # pylint: disable=too-many-statements
        ''' Check that filter patterns generate the expected python structures '''

        # simple =
        self.check_pattern("qual=4", 'qual = 4')
        # collapse == into =
        self.check_pattern("Qual == 50", 'qual = 50')

        # equivalent ORs
        self.check_pattern("qual >99 || qual <= -7", 'qual > 99 OR qual <= -7')
        self.check_pattern("qual >99 OR qual <= 7", 'qual > 99 OR qual <= 7')

        # equivalent ANDs
        self.check_pattern("qual >= 10 && qual < 50", 'qual >= 10 AND qual < 50')
        self.check_pattern("qual = 10 AND qual != 50", 'qual = 10 AND qual != 50')

        # valid function call
        self.check_pattern("MAX(7, 8)", 'MAX ( 7 , 8 )')

        # should fail, FOO() is not a function
        with self.assertRaises(ParseException):
            self.check_pattern("FOO(7, 8)", "FOO() is not a function")

        with self.assertRaises(ParseException):
            self.check_pattern("N_MISS > 17", "No missingness queries allowed")

        # FORMAT lookups. Note '' strings are generated.
        self.check_pattern("GT = '0/0'", "sample['GT'] = '0/0'")
        # smol fmt
        self.check_pattern('fmt/GT = "1/1"', "sample['GT'] = '1/1'")
        # BIG FMT, bare 0/1
        self.check_pattern("FMT/GT = 0/1", "sample['GT'] = '0/1'")
        # FORMAT w/ POS fixup
        self.check_pattern('FORMAT/GT="0/1", qual >= 73e-1 && pos > 8', "sample['GT'] = '0/1' OR qual >= 7.3 AND pos > 7")

        # FORMAT + INFO type lookups
        self.check_pattern("DP > 20", "CAST(sample['DP'] AS BIGINT) > 20")
        self.check_pattern("LAALTGC > 3e-2", "CAST(sample['LAALTGC'] AS DOUBLE) > 0.03")
        self.check_pattern("SVLEN > 50 || SVLEN < -50", "CAST(info['SVLEN'] AS BIGINT) > 50 OR CAST(info['SVLEN'] AS BIGINT) < -50")

        # mix FORMAT and INFO lookups
        self.check_pattern(""" GT='.|.' && SVTYPE = "INS" """, "sample['GT'] = '.|.' AND info['SVTYPE'] = 'INS'")
        self.check_pattern(""" DP > 20 && INFO/SVTYPE = "DEL" """, "CAST(sample['DP'] AS BIGINT) > 20 AND info['SVTYPE'] = 'DEL'")

        # other columns
        self.check_pattern("checkpoint = 5", "checkpoint = 5")
        self.check_pattern("reflen > 72", "reflen > 72")
        self.check_pattern("study_name = 'ajtrio'", "study_name = 'ajtrio'")
        self.check_pattern("aid = 'd177e9c3-1b86-4656-bba0-ea13386e1426'", "aid = 'd177e9c3-1b86-4656-bba0-ea13386e1426'")

        # VARID IS NULL
        self.check_pattern('ID = "."', "varid IS NULL")
        self.check_pattern("ID = '.'", "varid IS NULL")
        self.check_pattern("ID='.', FMT/GT == '0/1'", "varid IS NULL OR sample['GT'] = '0/1'")
        self.check_pattern("varid = '.'", "varid IS NULL")

        # These should all be sample['DP']
        self.check_pattern("DP = 10", "CAST(sample['DP'] AS BIGINT) = 10")
        self.check_pattern("((5 + 3) * 8) < DP", "( ( 5 + 3 ) * 8 ) < CAST(sample['DP'] AS BIGINT)")
        self.check_pattern("DP > ((5 + 3) * 8)", "CAST(sample['DP'] AS BIGINT) > ( ( 5 + 3 ) * 8 )")

        # genotype recognition
        self.check_pattern("GT != 1/1 and DP > 50", "sample['GT'] != '1/1' AND CAST(sample['DP'] AS BIGINT) > 50")
        self.check_pattern("GT != '1/1' and DP > 50", "sample['GT'] != '1/1' AND CAST(sample['DP'] AS BIGINT) > 50")
        self.check_pattern("GT == 0/1", "sample['GT'] = '0/1'")
        self.check_pattern("GT = ./.", "sample['GT'] = './.'")
        self.check_pattern("GT = 1|1", "sample['GT'] = '1|1'")

        # chrom fix helper
        self.check_pattern("chrom = '5'", "chrom = '5'")
        self.check_pattern('chrom = "5"', "chrom = '5'")
        self.check_pattern("ChRoM = 5", "chrom = '5'")
        self.check_pattern("CHROM = 'X'", "chrom = 'X'")
        self.check_pattern("CHROM = X", "chrom = 'X'")

        # PASS fix helper
        self.check_pattern("filter = 'PASS'", "filt = 'PASS'")
        self.check_pattern("FILTER = 'PASS'", "filt = 'PASS'")
        self.check_pattern("FILTER = PASS", "filt = 'PASS'")
        self.check_pattern("filt = lowq", "filt = 'lowq'")
        self.check_pattern("FILT = lowq", "filt = 'lowq'")

        # examples from bcftools. Uncommented examples are implemented, the rest are TODO.

        self.check_pattern("MIN(DV)>5", "MIN ( CAST(sample['DV'] AS BIGINT) ) > 5") # .. selects the whole site, evaluates min across all values and samples
        # SMPL_MIN(DV)>5  .. selects matching samples, evaluates within samples
        self.check_pattern("MIN(DV/DP)>0.3", "MIN ( CAST(sample['DV'] AS BIGINT) / CAST(sample['DP'] AS BIGINT) ) > 0.3")
        # self.check_pattern("MIN(DP)>10 & MIN(DV)>3", False) # TODO: support & |
        # FMT/DP>10  & FMT/GQ>10 .. both conditions must be satisfied within one sample
        self.check_pattern("FMT/DP>10 && FMT/GQ>10", "CAST(sample['DP'] AS BIGINT) > 10 AND CAST(sample['GQ'] AS BIGINT) > 10") # .. the conditions can be satisfied in different samples
        # QUAL>10 |  FMT/GQ>10   .. true for sites with QUAL>10 or a sample with GQ>10, but selects only samples with GQ>10
        self.check_pattern("QUAL>10 || FMT/GQ>10", "qual > 10 OR CAST(sample['GQ'] AS BIGINT) > 10") #   .. true for sites with QUAL>10 or a sample with GQ>10, plus selects all samples at such sites
        # TYPE="snp" && QUAL>=10 && (DP4[2]+DP4[3] > 2)
        # COUNT(GT="hom")=0      .. no homozygous genotypes at the site

        # This query: vvvvvvvvvv is correctly parsed but aggregates are not allowed in SELECTs in Athena
        self.check_pattern("AVG(GQ)>50", "AVG ( CAST(sample['GQ'] AS BIGINT) ) > 50") # .. average (arithmetic mean) of genotype qualities bigger than 50

        # ID=@file       .. selects lines with ID present in the file
        # ID!=@~/file    .. skip lines with ID present in the ~/file
        # MAF[0]<0.05    .. select rare variants at 5% cutoff

        # 0 vs. 1 based position translation
        self.check_pattern("POS>=100", "pos >= 99")
        self.check_pattern("pos > 10000 AND varend <= 22000", "pos > 9999 AND varend <= 22000")

    def test_missingness_patterns(self):
        ''' Check that filter patterns generate the expected python structures '''

        # SAMPLE_MISS
        self.check_pattern("SAMPLE_MISS < 0.1", 'S_MISS < 0.1', 'missingness')
        self.check_pattern("sample_miss >= 1", 'S_MISS >= 1', 'missingness')
        self.check_pattern("SAMPLE_MISSING < 0.3", 'S_MISS < 0.3', 'missingness')

        self.check_pattern("f_miss < 0.2", "CAST(infos['F_MISS'] AS DOUBLE) < 0.2", 'missingness')
        self.check_pattern("F_MISSING == 0.2", "CAST(infos['F_MISS'] AS DOUBLE) = 0.2", 'missingness')
        self.check_pattern("N_MISSING < 27", "CAST(infos['N_MISS'] AS DOUBLE) < 27", 'missingness')
        self.check_pattern("N_miss != 5", "CAST(infos['N_MISS'] AS DOUBLE) != 5", 'missingness')

        with self.assertRaises(ParseException):
            self.check_pattern("QUAL > 50", "not a missingness query", 'missingness')

        with self.assertRaises(ParseException):
            self.check_pattern("N_MISS > 3 OR F_MISS < 0.3", "one clause at a time", 'missingness')

        with self.assertRaises(ParseException):
            self.check_pattern("(SAMPLE_MISS > 0.7)", "no parens", 'missingness')

if __name__ == '__main__':
    unittest.main(verbosity=2)
