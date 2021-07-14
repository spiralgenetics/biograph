# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph


class DnaTestCases(unittest.TestCase):

    def test_basic(self):
        seq = biograph.Sequence("GATTACA")
        self.assertEqual(seq.rev_comp(), biograph.Sequence("TGTAATC"))
        self.assertGreater(seq, biograph.Sequence("GATTAC"))
        self.assertLess(seq, biograph.Sequence("GATTACC"))
        self.assertEqual(len(seq), 7)
        self.assertEqual(seq[::], "GATTACA")
        self.assertEqual(seq[0], biograph.Sequence("G"))
        self.assertEqual(seq[0], "G")
        self.assertEqual(seq[:3], biograph.Sequence("GAT"))
        self.assertEqual(seq[:3], "GAT")
        self.assertEqual(seq[:3:1], "GAT")
        self.assertEqual(seq[-1], "A")
        self.assertEqual(seq[-len(seq)], "G")
        self.assertEqual(str(seq), "GATTACA")
        self.assertEqual(seq, eval(repr(seq))) # pylint: disable=eval-used

    def test_bad_range(self):
        seq = biograph.Sequence("GATTACA")
        with self.assertRaises(IndexError):
            _ = seq[::-1]
        with self.assertRaises(IndexError):
            _ = seq[::2]
        with self.assertRaises(IndexError):
            _ = seq[len(seq)]
        with self.assertRaises(IndexError):
            _ = seq[-(len(seq) + 1)]

    def test_modify(self):
        seq = biograph.Sequence("GATTACA")
        seq[3] = "A"
        self.assertEqual(seq, "GATAACA")
        seq += "T"
        self.assertEqual(seq, "GATAACAT")
        seq = seq + "A"
        self.assertEqual(seq, "GATAACATA")
        seq += biograph.Sequence("T")
        self.assertEqual(seq, "GATAACATAT")
        seq = seq + biograph.Sequence("A")
        self.assertEqual(seq, "GATAACATATA")
        seq = biograph.Sequence("C") + seq
        self.assertEqual(seq, "CGATAACATATA")
        seq = "G" + seq
        self.assertEqual(seq, "GCGATAACATATA")
        seq[-3] = "C"
        self.assertEqual(seq, "GCGATAACATCTA")

    def test_modify_range(self):
        for start, end, repl in [
                (3, 5, "A"),
                (3, -1, "T"),
                (-5, -3, "C")]:
            print("Attempting to replace %d : %d = %s" % (start, end, repl))
            seq = biograph.Sequence("GATTACA")
            expected = list("GATTACA")
            print("Before replace: %s" % seq)
            seq[start:end] = repl
            print("After replace: %s" % seq)
            expected[start:end] = list(repl)
            self.assertEqual(seq, "".join(expected))
        seq = biograph.Sequence("GATTACA")
        seq[::] = "CAT"
        self.assertEqual(seq, "CAT")

    def test_bad_modify_range(self):
        seq = biograph.Sequence("GATTACA")
        with self.assertRaises(IndexError):
            _ = seq[::-1] = "A"
        with self.assertRaises(IndexError):
            _ = seq[::2] = "A"
        with self.assertRaises(IndexError):
            _ = seq[len(seq)] = "A"
        with self.assertRaises(IndexError):
            _ = seq[-(len(seq) + 1)] = "A"

    def test_iterate(self):
        seq = biograph.Sequence("GATTACA")
        result = biograph.Sequence()
        for base in seq:
            self.assertEqual(1, len(base))
            result = result + base
        self.assertEqual(seq, result)

if __name__ == '__main__':
    unittest.main(verbosity=2)
