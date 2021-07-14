# pylint: disable=missing-docstring

import unittest
import gzip
import tempfile

from biograph.vdb import vcf_to_parquet, anno_to_parquet
import biograph.vdb.parquet
import pyarrow as pa
import pandas as pd

class ParquetTestCases(unittest.TestCase):
    def setUp(self):
        # Make sure we test having multiple input chunks
        biograph.vdb.parquet.ConverterBase.CHUNK_ROWS = 50
        self.maxDiff = None
        self.span_bases = 59

    def get_parquet_rows(self, pq_path):
        "Returns all records from the given parquet file as a list of dict objects"
        f = pa.parquet.ParquetDataset(pq_path)
        # Order within individual .parquet files is sorted by pos, but order is not
        # guaranteed for a whole partitioned dataset, so sort here.
        df = f.read().to_pandas().sort_values(["chrom", "pos"])
        # Convert NAN to Python's 'None':
        df = df.where(pd.notnull(df), None)
        records = df.to_dict(orient='records')

        # Check proper 'spans' value remove it.
        for row in records:
            self.assertLessEqual(row['pos'], row['varend'])

            expected_spans = []
            span = row['pos'] // self.span_bases
            while span <= row['varend'] // self.span_bases:
                expected_spans.append(span)
                span += 1
            self.assertListEqual(expected_spans, list(row['spans']), msg=repr(row))
            del row['spans']

        return records

    # @unittest.skip(True)
    def test_vcf(self):
        with tempfile.TemporaryDirectory() as td:
            out_fn = td
            with gzip.open("golden/ftest/vdb/vdb002.vcf.gz", "rb") as fh:
                header_lines = vcf_to_parquet("grch37", fh, out_fn, span_bases=self.span_bases)
            self.assertEqual(len(header_lines), 140)
            self.assertEqual(header_lines[0], '##fileformat=VCFv4.1')
            actual = self.get_parquet_rows(out_fn)
            self.assertEqual(len(actual), 2500)
            first_row = actual[0]
            self.assertDictEqual(first_row, {
                'alt': 'A',
                'chrom': '1',
                'filt': 'lowq',
                'info': [('NS', '1'), ('POP', ''), ('SVLEN', '-37'), ('SVTYPE', 'DEL')],
                'pos': 10402,
                'qual': 11.0,
                'ref': 'ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC',
                'reflen': 38,
                'sample': [('GT', '0/1'),
                           ('AC', '6,8'),
                           ('AD', '7,1'),
                           ('DC', '6,8'),
                           ('DCC', '1,0'),
                           ('DDC', '5,8'),
                           ('DMO', '133,143'),
                           ('DP', '8'),
                           ('DS', '45,18'),
                           ('DXO', '148,148'),
                           ('EC', '6,8'),
                           ('GQ', '1'),
                           ('LAALTGC', '0.68'),
                           ('LAALTSEQLEN', '125'),
                           ('LALANCH', '0'),
                           ('LARANCH', '0'),
                           ('LAREFGC', '0.641975'),
                           ('LAREFSPAN', '162'),
                           ('LASCORE', '0'),
                           ('MC', '6,8'),
                           ('MO', '133,143'),
                           ('MP', '0,0'),
                           ('NR', '8,8'),
                           ('NUMASM', '1'),
                           ('OV', '0'),
                           ('PAD', '0,0'),
                           ('PDP', '0'),
                           ('PG', '0|1'),
                           ('PI', '7854514'),
                           ('PL', '1,3,16'),
                           ('RC', '0'),
                           ('UC', '2,8'),
                           ('UCC', '0,0'),
                           ('UDC', '2,8'),
                           ('UMO', '146,143'),
                           ('US', '10,18'),
                           ('UXO', '146,148'),
                           ('XC', '6,8'),
                           ('XO', '148,148')],
                'varend': 10440,
                'varid': None
            })

    # @unittest.skip(True)
    def test_vcf_verbose_threads(self):
        with tempfile.TemporaryDirectory() as td:
            out_fn = td
            with gzip.open("golden/ftest/vdb/vdb002.vcf.gz", "rb") as fh:
                header_lines = vcf_to_parquet("grch37", fh, out_fn, nthreads=10, span_bases=self.span_bases, verbose=True)
            self.assertEqual(len(header_lines), 140)
            self.assertEqual(header_lines[0], '##fileformat=VCFv4.1')
            actual = self.get_parquet_rows(out_fn)
            self.assertEqual(len(actual), 2500)
            first_row = actual[0]
            self.assertDictEqual(first_row, {
                'alt': 'A',
                'chrom': '1',
                'filt': 'lowq',
                'info': [('NS', '1'), ('POP', ''), ('SVLEN', '-37'), ('SVTYPE', 'DEL')],
                'pos': 10402,
                'qual': 11.0,
                'ref': 'ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC',
                'reflen': 38,
                'sample': [('GT', '0/1'),
                           ('AC', '6,8'),
                           ('AD', '7,1'),
                           ('DC', '6,8'),
                           ('DCC', '1,0'),
                           ('DDC', '5,8'),
                           ('DMO', '133,143'),
                           ('DP', '8'),
                           ('DS', '45,18'),
                           ('DXO', '148,148'),
                           ('EC', '6,8'),
                           ('GQ', '1'),
                           ('LAALTGC', '0.68'),
                           ('LAALTSEQLEN', '125'),
                           ('LALANCH', '0'),
                           ('LARANCH', '0'),
                           ('LAREFGC', '0.641975'),
                           ('LAREFSPAN', '162'),
                           ('LASCORE', '0'),
                           ('MC', '6,8'),
                           ('MO', '133,143'),
                           ('MP', '0,0'),
                           ('NR', '8,8'),
                           ('NUMASM', '1'),
                           ('OV', '0'),
                           ('PAD', '0,0'),
                           ('PDP', '0'),
                           ('PG', '0|1'),
                           ('PI', '7854514'),
                           ('PL', '1,3,16'),
                           ('RC', '0'),
                           ('UC', '2,8'),
                           ('UCC', '0,0'),
                           ('UDC', '2,8'),
                           ('UMO', '146,143'),
                           ('US', '10,18'),
                           ('UXO', '146,148'),
                           ('XC', '6,8'),
                           ('XO', '148,148')],
                'varend': 10440,
                'varid': None
            })

    # @unittest.skip(True)
    def parse_anno(self, fmt, build, filename):
        with tempfile.TemporaryDirectory() as td:
            out_fn = td
            with open(filename, "rb") as fh:
                header_lines = anno_to_parquet(fmt, build, fh, out_fn, span_bases=self.span_bases)
            records = self.get_parquet_rows(out_fn)
            return header_lines, records

    # @unittest.skip(True)
    def test_vcf_anno(self):
        header_lines, actual = self.parse_anno("vcf", "grch37", "golden/ftest/vdb/dbsnp-sample.vcf")
        self.assertEqual(len(header_lines), 85)
        self.assertEqual(header_lines[0], '##fileformat=VCFv4.0')
        self.assertEqual(len(actual), 4162)
        first_row = actual[0]
        self.assertDictEqual(first_row,
                             {'reflen': 1, 'chrom': 'X', 'pos': 60001, 'varend': 60002,
                              'varid': 'rs1226858834', 'ref': 'T', 'alt': 'A', 'qual': None,
                              'filt': None,
                              'info': [('RS', '1226858834'), ('RSPOS', '60002'), ('dbSNPBuildID', '151'), ('SSR', '0'), ('SAO', '0'), ('VP', '0x050000000005000002000100'), ('WGT', '1'), ('VC', 'SNV'), ('ASP', ''), ('TOPMED', '0.99631275484199796,0.00368724515800203')],
                              'source': None, 'feature': None, 'score': None,
                              'frame': None, 'strand': None, 'attributes': None})

    # @unittest.skip(True)
    def test_gtf_anno(self): # GTF, also known as GFFv2
        header_lines, actual = self.parse_anno("gtf", "grch37", "golden/ftest/vdb/gtf-37-sample.gtf")
        self.assertEqual(len(header_lines), 5)
        self.assertEqual(header_lines[0], '#!genome-build GRCh37.p13')
        self.assertEqual(len(actual), 7)
        first_gene_row = [row for row in actual if row['feature'] == 'gene'][0]
        self.assertDictEqual(first_gene_row,
                             {'reflen': 2543, 'chrom': '1', 'pos': 11868, 'varend': 14411,
                              'varid': 'DDX11L1', 'ref': None, 'alt': None, 'qual': None,
                              'filt': None, 'info': None, 'source': 'pseudogene',
                              'feature': 'gene', 'score': None, 'frame': None, 'strand': '+',
                              'attributes': [('gene_id', 'ENSG00000223972'),
                                             ('gene_name', 'DDX11L1'),
                                             ('gene_source', 'ensembl_havana'),
                                             ('gene_biotype', 'pseudogene')]})

    # @unittest.skip(True)
    def test_gff_anno(self): # GFFv3
        header_lines, actual = self.parse_anno("gff", "grch38", "golden/ftest/vdb/gff-38-sample.gff3")
        self.assertEqual(len(header_lines), 200)
        self.assertEqual(header_lines[0], '##gff-version 3')
        self.assertEqual(len(actual), 64)

        third_row = actual[2]
        # Floats are only almost equal.
        self.assertAlmostEqual(third_row['score'], 0.9990000128746033)
        del third_row['score']

        self.assertDictEqual(third_row, {
            'alt': None,
            'attributes': [('logic_name', 'eponine')],
            'chrom': '1',
            'feature': 'biological_region',
            'filt': None,
            'frame': None,
            'info': None,
            'pos': 10649,
            'qual': None,
            'ref': None,
            'reflen': 7,
            'source': '.',
            'strand': '+',
            'varend': 10656,
            'varid': 'eponine'
        })

        other_rows = [a for a in actual if a['varid'] == 'DDX11L1']
        self.assertEqual(len(other_rows), 1)
        self.assertDictEqual(other_rows[0],
                             {
                                 'alt': None,
                                 'attributes': [('ID', 'gene:ENSG00000223972'),
                                                ('Name', 'DDX11L1'),
                                                ('biotype', 'transcribed_unprocessed_pseudogene'),
                                                ('description',
                                                 'DEAD/H-box helicase 11 like 1 (pseudogene) [Source:HGNC Symbol;Acc:HGNC:37102]'),
                                                ('gene_id', 'ENSG00000223972'),
                                                ('logic_name', 'havana_homo_sapiens'),
                                                ('version', '5')],
                                 'chrom': '1',
                                 'feature': 'pseudogene',
                                 'filt': None,
                                 'frame': None,
                                 'info': None,
                                 'pos': 11868,
                                 'qual': None,
                                 'ref': None,
                                 'reflen': 2540,
                                 'score': None,
                                 'source': 'havana',
                                 'strand': '+',
                                 'varend': 14408,
                                 'varid': 'DDX11L1'
                             })

if __name__ == '__main__':
    unittest.main(verbosity=2)
