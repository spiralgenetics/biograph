"""
refhash: Compute contig:len hashes

//python/functest:refhash_test

"""
import unittest

from python.functest.utils.setup import (
    ftest_setup,
    ftest_teardown,
    ftest_module_setup
)
from python.functest.utils.defaults import (
    GOLDEN_DIR
)
from biograph.tools.refhash import refhash, refstyle
from biograph.utils import cmd_exe

def setUpModule():
    """ Announce ourselves by printing the module docstring. """
    print(__doc__)

    # Module requirements
    ftest_module_setup()

class RefHashTestCases(unittest.TestCase):
    """Test formats"""

    data_dir = None

    def setUp(self):
        ftest_setup(self)

    def tearDown(self):
        ftest_teardown(self)

    def test_refhash_vcf(self):
        """Test VCF """
        # vcf
        rh = refhash(f'{GOLDEN_DIR}/e_coli_sas.vcf')
        self.assertEqual(rh.digest(), '4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432')
        self.assertEqual(rh.common_name(), 'e_coli_k12_ASM584v1')
        self.assertEqual(rh.verbose_name(), 'refhash=4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432,name=e_coli_k12_ASM584v1')

        # gzvcf
        rh = refhash(f'{GOLDEN_DIR}/vdb/vdb002.vcf.gz')
        self.assertEqual(rh.digest(), '1e4ef0c15393ae133ad336a36a376bb62e564a43a65892966002e75713282aec')
        self.assertEqual(rh.common_name(), 'hs37d5')
        self.assertEqual(rh.verbose_name(), 'refhash=1e4ef0c15393ae133ad336a36a376bb62e564a43a65892966002e75713282aec,name=hs37d5')

    def test_refhash_fasta(self):
        """Test FASTA """
        # e.coli
        rh = refhash(f'datasets/reference/e_coli_k12_ASM584v1/source.fasta')
        self.assertEqual(rh.digest(), '4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432')
        self.assertEqual(rh.common_name(), 'e_coli_k12_ASM584v1')
        self.assertEqual(rh.verbose_name(), 'refhash=4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432,name=e_coli_k12_ASM584v1')

        # unknown
        rh = refhash(f'{GOLDEN_DIR}/E_coli_IS_seqs.fasta.gz')
        self.assertEqual(rh.digest(), 'bb16ab78764511b821e36055f4582f306b4c362ec2f6fef9c8dabd392c081b5e')
        self.assertEqual(rh.common_name(), 'unknown-bb16ab78')
        self.assertEqual(rh.verbose_name(), 'refhash=bb16ab78764511b821e36055f4582f306b4c362ec2f6fef9c8dabd392c081b5e,name=unknown-bb16ab78')

        # some lines have LN:
        rh = refhash(f'{GOLDEN_DIR}/lntest.fasta')
        self.assertEqual(rh.digest(), 'c8fe760a958bbaeba1a602475bedead8c66d404462f2c8f6c7ffb46bc6347b66')
        self.assertEqual(rh.common_name(), 'unknown-c8fe760a')
        self.assertEqual(rh.verbose_name(), 'refhash=c8fe760a958bbaeba1a602475bedead8c66d404462f2c8f6c7ffb46bc6347b66,name=unknown-c8fe760a')

    def test_refhash_bg_refdir(self):
        """Test bg refdir """
        # e.coli (refhash should match fasta test)
        rh = refhash(f'datasets/reference/e_coli_k12_ASM584v1/')
        self.assertEqual(rh.digest(), '4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432')
        self.assertEqual(rh.common_name(), 'e_coli_k12_ASM584v1')
        self.assertEqual(rh.verbose_name(), 'refhash=4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432,name=e_coli_k12_ASM584v1')

        # refdir
        rh = refhash(f'/reference/grch38/')
        self.assertEqual(rh.digest(), '1f5faf40c2b1b8715e9df75375cb392117a9c5734fca790e6399d7a50e90ebdd')
        self.assertEqual(rh.common_name(), 'grch38')
        self.assertEqual(rh.verbose_name(), 'refhash=1f5faf40c2b1b8715e9df75375cb392117a9c5734fca790e6399d7a50e90ebdd,name=grch38')

        # karyotype.json
        rh = refhash(f'/reference/hs37d5/karyotype.json')
        self.assertEqual(rh.digest(), '1e4ef0c15393ae133ad336a36a376bb62e564a43a65892966002e75713282aec')
        self.assertEqual(rh.common_name(), 'hs37d5')
        self.assertEqual(rh.verbose_name(), 'refhash=1e4ef0c15393ae133ad336a36a376bb62e564a43a65892966002e75713282aec,name=hs37d5')

    def test_refhash_bam_header(self):
        """Test sam/bam/cram header """
        header = cmd_exe("samtools view -H datasets/bams/e_coli/e_coli_test.bam").stdout

        with open("header.txt", "w") as f:
            f.write(header)

        rh = refhash(f'header.txt')
        self.assertEqual(rh.digest(), '4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432')
        self.assertEqual(rh.common_name(), 'e_coli_k12_ASM584v1')
        self.assertEqual(rh.verbose_name(), 'refhash=4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432,name=e_coli_k12_ASM584v1')

    def test_refhash_fail(self):
        """Test failure modes """

        # no_contigs_present
        with self.assertRaises(SystemExit):
            refhash(f'/dev/null')

        # dir, not bgdir
        with self.assertRaises(SystemExit):
            refhash(f'/tmp/')

        # no contigs in vcf
        rh = refhash(f'{GOLDEN_DIR}/microcontigs.vcf.gz')
        self.assertEqual(rh.verbose_name(), 'refhash=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855,name=no_contigs_present')

        # no file
        with self.assertRaises(FileNotFoundError):
            refhash(f'/file/not/found')

        # invalid fasta
        with open('invalid.fa', 'w') as f:
            f.write('>')
        with self.assertRaises(SystemExit):
            refhash(f'invalid.fa')

    def test_name_translation(self):
        """ to_ebi() etc. """
        chroms = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
                  '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
                  'X', 'Y', 'MT')

        rh = refhash()

        for build in ('GRCh37', 'GRCh38'):
            # Invalid to_ebi lookups are simply returned
            self.assertEqual('invalid', rh.to_ebi('invalid', build))

            for chrom in (chroms):
                self.assertEqual(chrom, rh.to_ebi(rh.to_genbank(chrom, build), build))
                self.assertEqual(chrom, rh.to_ebi(rh.to_refseq(chrom, build), build))
                self.assertEqual(chrom, rh.to_ebi(rh.to_ucsc(chrom, build), build))

            for style in (refstyle.UCSC, refstyle.GenBank, refstyle.RefSeq):
                # Invalid from_ebi lookups are simply returned
                self.assertEqual('invalid', rh.from_ebi('invalid', build, style))

                for chrom in (chroms):
                    new_chrom = rh.from_ebi(chrom, build, style)
                    self.assertNotEqual(chrom, new_chrom)
                    self.assertEqual(chrom, rh.to_ebi(new_chrom, build))

                self.assertEqual(rh.to_native(chrom, build, style), rh.from_ebi(chrom, build, style))

if __name__ == '__main__':
    unittest.main(verbosity=2)
