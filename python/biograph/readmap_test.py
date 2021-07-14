# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph


class ReadmapTestCases(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("golden/e_coli_merged.bg")
        cls.seqset = cls.bg.seqset
        cls.rm1 = cls.bg.open_readmap("test_accession_id")
        cls.rm2 = cls.bg.open_readmap("e_coli_test")

        cls.unmerged_bg = biograph.BioGraph("golden/e_coli_10000snp.bg")
        cls.unmerged_seqset = cls.unmerged_bg.seqset
        cls.unmerged_rm = cls.unmerged_bg.open_readmap("")

    def test_size(self):
        self.assertEqual(self.rm1.size(), 16888)
        self.assertEqual(self.rm2.size(), 76094)
        self.assertEqual(self.unmerged_rm.size(), 16888)

    def test_read_lens(self):
        self.assertEqual(self.rm1.min_read_len(), 30)
        self.assertEqual(self.rm2.min_read_len(), 100)
        self.assertEqual(self.rm1.max_read_len(), 35)
        self.assertEqual(self.rm2.max_read_len(), 100)

    def test_read_get_by_id(self):
        read = self.rm1.get_read_by_id(123)
        self.assertEqual(read.get_read_id(), 123)
        self.assertEqual(len(read), 35)
        self.assertEqual(read.get_seqset_entry().sequence(),
                         biograph.Sequence("AAAATTACAGAGTACACAACATCCATGAAACGCAT"))

    def test_get_empty_prefixes(self):
        reads = self.rm1.get_prefix_reads(self.seqset.find("A"))
        read_count = 0
        for _ in reads:
            read_count = read_count + 1
        self.assertEqual(read_count, 0)

    def test_read_get_by_id_out_of_range(self):
        self.rm1.get_read_by_id(16887)
        with self.assertRaises(RuntimeError):
            self.rm1.get_read_by_id(16888)

    def test_get_prefixes_equal(self):
        r = self.seqset.find("AAAATTACAGAGTACACAACATCCATGAAACGCAT")
        self.assertFalse(list(self.rm2.get_prefix_reads(r)))
        reads = list(self.rm1.get_prefix_reads(r))
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].get_read_id(), 123)
        self.assertTrue(reads[0].is_original_orientation())
        self.assertEqual(reads[1].get_read_id(), 124)
        self.assertFalse(reads[1].is_original_orientation())
        for read in reads:
            self.assertEqual(read.get_seqset_entry().sequence(), r.sequence())
        print(self.rm2)

    def test_get_prefixes_prefix(self):
        r = self.seqset.find("GGATGAAATGAGTTGCCATCTGGTGCTCACCACGG")
        self.assertFalse(list(self.rm2.get_prefix_reads(r)))
        print(self.rm2)
        print(self.rm1)
        print("Starting for")
        for x in self.rm1.get_prefix_reads(r):
            print("Printing for")
            print(x)
            print("Done printing for")
        print("Done for")
        reads = [(read.get_read_id(), str(read.get_seqset_entry().sequence()))
                 for read in self.rm1.get_prefix_reads(r)]
        print(self.rm2)
        print(self.rm1)
        self.maxDiff = 100000
        self.assertEqual(reads,
                         [(11160, "GGATGAAATGAGTTGCCATCTGGTGCTCACCACGG"),
                          (11161, "GGATGAAATGAGTTGCCATCTGGTGCTCAC"),
                          (11162, "GGATGAAATGAGTTGCCATCTGGTGCTCACC"),
                          (11163, "GGATGAAATGAGTTGCCATCTGGTGCTCACCA"),
                          (11164, "GGATGAAATGAGTTGCCATCTGGTGCTCACCACG"),
                          (11165, "GGATGAAATGAGTTGCCATCTGGTGCTCACCACGG")])
        print(self.rm2)
        print(self.rm1)

    def test_get_mate(self):
        # 63816 is part of a pair
        read = self.rm2.get_read_by_id(63816)
        self.assertTrue(read.has_mate())
        self.assertEqual(read.get_mate().get_read_id(), 37698)

        # 35906 is not part of a pair
        read = self.rm2.get_read_by_id(35906)
        self.assertFalse(read.has_mate())
        with self.assertRaises(RuntimeError):
            read.get_mate()

    def test_wrong_seqset(self):
        r = self.unmerged_seqset.find("AAAATTACAGAGTACACAACATCCATGAAACGCAT")
        with self.assertRaises(RuntimeError):
            self.rm1.get_prefix_reads(r)

    def test_reads_containing(self):
        r = self.seqset.find("TGAAATGAGTTGCCATCTGGTGCTCAC")
        reads = [(offset, read.get_read_id(), str(read.get_seqset_entry().sequence()))
                 for offset, read in self.rm1.get_reads_containing(r)]
        self.assertEqual(
            reads, [(0, 14544, 'TGAAATGAGTTGCCATCTGGTGCTCACCACGG'),
                    (1, 3553, 'ATGAAATGAGTTGCCATCTGGTGCTCACCACGG'),
                    (3, 11160, 'GGATGAAATGAGTTGCCATCTGGTGCTCACCACGG'),
                    (3, 11161, 'GGATGAAATGAGTTGCCATCTGGTGCTCAC'),
                    (3, 11162, 'GGATGAAATGAGTTGCCATCTGGTGCTCACC'),
                    (3, 11163, 'GGATGAAATGAGTTGCCATCTGGTGCTCACCA'),
                    (3, 11164, 'GGATGAAATGAGTTGCCATCTGGTGCTCACCACG'),
                    (3, 11165, 'GGATGAAATGAGTTGCCATCTGGTGCTCACCACGG')]
        )

    def test_pair_stats(self):
        st1 = self.rm1.get_pair_stats()
        self.assertEqual(st1.paired_reads, 0)
        self.assertEqual(st1.paired_bases, 0)
        self.assertEqual(st1.unpaired_reads, 8444)
        self.assertEqual(st1.unpaired_bases, 288464)

        st2 = self.rm2.get_pair_stats()
        self.assertEqual(st2.paired_reads, 25762)
        self.assertEqual(st2.paired_bases, 2576200)
        self.assertEqual(st2.unpaired_reads, 12285)
        self.assertEqual(st2.unpaired_bases, 1228500)

if __name__ == '__main__':
    unittest.main(verbosity=2)
