# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph


class SeqsetTestCases(unittest.TestCase):

    def setUp(self):
        self.bg = biograph.BioGraph("golden/e_coli_10000snp.bg")
        self.seqset = self.bg.seqset
        self.gattaca = self.seqset.find(biograph.Sequence("GATTACA"))

    def test_size(self):
        self.assertEqual(self.seqset.size(), 19935)

    def test_empty_entry(self):
        empty = self.seqset.empty_entry()
        self.assertTrue(empty is not None)
        self.assertEqual(len(empty), 0)
        self.assertEqual(empty.sequence(), biograph.Sequence(""))
        self.assertEqual(str(empty), "<SeqsetEntry 0-19935: >")
        self.assertEqual(repr(empty), "<SeqsetEntry 0-19935: >")
        self.assertEqual(empty, self.seqset.find(biograph.Sequence("")))
        self.assertEqual(empty.get_begin_entry_id(), 0)
        self.assertEqual(empty.get_end_entry_id(), self.seqset.size())

    def test_find_existing(self):
        self.assertTrue(self.gattaca is not None)
        self.assertEqual(len(self.gattaca), len("GATTACA"))
        self.assertEqual(
            self.gattaca.sequence(), biograph.Sequence("GATTACA"))
        self.assertEqual(self.gattaca.sequence(3), biograph.Sequence("GAT"))
        self.assertEqual(
            str(self.gattaca), "<SeqsetEntry 11076-11077: GATTACA>")
        self.assertEqual(
            repr(self.gattaca), "<SeqsetEntry 11076-11077: GATTACA>")
        self.assertEqual(self.gattaca.get_begin_entry_id(), 11076)
        self.assertEqual(self.gattaca.get_end_entry_id(), 11077)

    def test_find_not_existing(self):
        not_found = self.seqset.find(biograph.Sequence("GATTACAGATTACA"))
        self.assertTrue(not_found is None)

    def test_truncate(self):
        truncated = self.gattaca.truncate(8)
        self.assertTrue(truncated == self.gattaca)
        self.assertFalse(truncated != self.gattaca)
        truncated = self.gattaca.truncate(7)
        self.assertTrue(truncated == self.gattaca)
        self.assertFalse(truncated != self.gattaca)
        truncated = self.gattaca.truncate(6)
        self.assertFalse(truncated == self.gattaca)
        self.assertTrue(truncated != self.gattaca)
        self.assertEqual(str(truncated.sequence()), "GATTAC")

    def test_sequence(self):
        self.assertEqual(str(self.gattaca.sequence(0)), "")
        self.assertEqual(str(self.gattaca.sequence(1)), "G")
        self.assertEqual(str(self.gattaca.sequence(2)), "GA")
        self.assertEqual(str(self.gattaca.sequence(3)), "GAT")
        self.assertEqual(str(self.gattaca.sequence(4)), "GATT")
        self.assertEqual(str(self.gattaca.sequence(5)), "GATTA")
        self.assertEqual(str(self.gattaca.sequence(6)), "GATTAC")
        self.assertEqual(str(self.gattaca.sequence(7)), "GATTACA")
        self.assertEqual(str(self.gattaca.sequence(8)), "GATTACA")

    def test_pop_front(self):
        self.assertEqual(str(self.gattaca.pop_front().sequence()), "ATTACA")
        r = self.seqset.find(biograph.Sequence("T")).pop_front()
        self.assertEqual(r, self.seqset.empty_entry())
        with self.assertRaises(RuntimeError):
            r.pop_front()

    def test_get_entry_by_id(self):
        with self.assertRaises(TypeError):
            self.seqset.get_entry_by_id(-1)
        by_id = self.seqset.get_entry_by_id(0)
        self.assertEqual(by_id.sequence(), biograph.Sequence("AAAAAAAAGCCCGCACTGTCAGGTGCGGGC"))
        by_id = self.seqset.get_entry_by_id(6645)
        self.assertEqual(by_id.sequence(), biograph.Sequence("CCCCACCTGACGCGAAAGGATCGCAGT"))
        self.assertEqual(by_id.get_begin_entry_id(), 6645)
        self.assertEqual(by_id.get_end_entry_id(), 6646)
        by_id = self.seqset.get_entry_by_id(19934)
        self.assertEqual(by_id.sequence(), biograph.Sequence("TTTTTTTTCGACCAAAGGTAACGAGGTAACAACCA"))
        with self.assertRaises(RuntimeError):
            self.seqset.get_entry_by_id(19935)

    def test_push_front(self):
        pushed = self.gattaca.push_front("G")
        self.assertTrue(pushed is None)
        pushed = self.gattaca.push_front("T")
        self.assertEqual(pushed, self.seqset.find("TGATTACA"))

    def test_push_front_drop(self):
        pushed = self.gattaca.push_front_drop("G")
        self.assertEqual(pushed.sequence(), biograph.Sequence("GGATTAC"))
        pushed = self.gattaca.push_front_drop("G", 8)
        self.assertTrue(pushed is None)
        pushed = self.gattaca.push_front_drop("G", 7)
        self.assertTrue(pushed is None)
        pushed = self.gattaca.push_front_drop("G", 6)
        self.assertEqual(pushed.sequence(), biograph.Sequence("GGATTAC"))

if __name__ == '__main__':
    unittest.main(verbosity=2)
