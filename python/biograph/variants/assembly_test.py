# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import itertools

import biograph
import biograph.variants as bgexvar

class AssemblyTestCases(unittest.TestCase):
    def test_asm_construct(self):
        asm = bgexvar.Assembly(1000, 1050, biograph.Sequence("GATTACA"), 5)
        self.assertEqual(asm.seq, "GATTACA")

    def test_asm_construct_autoconvert(self):
        asm = bgexvar.Assembly(1000, 1050, "GATTACA", 5)
        self.assertEqual(asm.seq, "GATTACA")

    def test_asm_fields(self): # pylint:disable=too-many-statements
        asm = bgexvar.Assembly(1000, 1050, biograph.Sequence("GATTACA"), 5)
        self.assertEqual(asm.seq, "GATTACA")
        self.assertEqual(asm.assembly_id, 5)
        self.assertEqual(asm.left_offset, 1000)
        self.assertEqual(asm.right_offset, 1050)
        self.assertIsNone(asm.edge_coverage)
        self.assertIsNone(asm.read_coverage)
        self.assertEqual(asm.matches_reference, False)
        self.assertEqual(asm.min_overlap, 0)
        self.assertEqual(asm.generated_by, "")
        self.assertEqual(asm.tags, [])

        asm.seq = "ATTAC"
        self.assertEqual(asm.seq, "ATTAC")
        asm.left_offset = 1234
        self.assertEqual(asm.left_offset, 1234)
        asm.right_offset = 5678
        self.assertEqual(asm.right_offset, 5678)

        asm.assembly_id = 9012
        self.assertEqual(asm.assembly_id, 9012)

        asm.min_overlap = 5
        self.assertEqual(asm.min_overlap, 5)
        asm.generated_by = "PYTHON"
        self.assertEqual(asm.generated_by, "PYTHON")
        self.assertEqual(asm.tags, ["PYTHON"])

        asm.tags = ["A", "B"]
        self.assertEqual(asm.generated_by, "A,B")
        self.assertEqual(asm.tags, ["A", "B"])

        ReadCoverageRead = bgexvar.ReadCoverageRead
        self.assertEqual(repr(ReadCoverageRead(1, 2, 3)), "ReadCoverageRead(1, 2, 3)")

        asm.read_coverage = bgexvar.ReadCoverage(len("GATTACA"), [])
        self.assertEqual(asm.read_coverage, bgexvar.ReadCoverage(len("GATTACA"), []))
        cov_table = asm.read_coverage
        self.assertEqual(len(cov_table), 0)
        self.assertEqual(list(cov_table), [])
        self.assertEqual(cov_table.assembly_len(), len("GATTACA"))
        readlist = list(cov_table)
        readlist.append(ReadCoverageRead(1, 2, 3))
        asm.read_coverage = bgexvar.ReadCoverage(len("GATTACA"), readlist)
        cov_table = asm.read_coverage
        self.assertEqual(len(cov_table), 1)
        self.assertEqual(list(cov_table), [ReadCoverageRead(1, 2, 3)])
        readlist = list(cov_table)
        readlist.append(ReadCoverageRead(1, 2, 3))
        asm.read_coverage = bgexvar.ReadCoverage(len("GATTACA"), readlist)
        cov_table = asm.read_coverage
        self.assertEqual(list(cov_table), [ReadCoverageRead(1, 2, 3)])
        self.assertEqual(cov_table, bgexvar.ReadCoverage(len("GATTACA"), [ReadCoverageRead(1, 2, 3)]))
        self.assertNotEqual(cov_table, bgexvar.ReadCoverage(12345, [ReadCoverageRead(1, 2, 3)]))
        self.assertNotEqual(cov_table, bgexvar.ReadCoverage(len("GATTACA"), [ReadCoverageRead(2, 1, 5)]))
        # These are not convertable to ReadCoverageReads.
        with self.assertRaises(TypeError):
            bgexvar.ReadCoverage([1, 2, 3])
        # Make sure generator exceptions are propagated.
        with self.assertRaises(KeyError):
            bgexvar.ReadCoverage(len("GATTACA"),
                                 ({0: ReadCoverageRead(1, 2, 3), 2: ReadCoverageRead(2, 1, 5)}[x]
                                  for x in range(3)))

    def test_asm_phase_ids(self):
        asm = bgexvar.Assembly(1000, 1050, biograph.Sequence("GATTACA"), 5)
        self.assertSequenceEqual(asm.phase_ids, [])
        self.assertFalse(asm.phase_ids)
        asm.phase_ids = bgexvar.PhaseSet(["c", "a", "b", "c"])
        # phase_ids becomes sorted and deduplicated:
        self.assertSequenceEqual(list(asm.phase_ids), ["a", "b", "c"])
        self.assertTrue(asm.phase_ids)
        self.assertEqual(str(asm.phase_ids), "(a,b,c)")
        self.assertEqual(repr(asm.phase_ids), "PhaseSet([\"a\",\"b\",\"c\"])")
        asm.phase_ids.clear()
        self.assertSequenceEqual(list(asm.phase_ids), [])
        self.assertFalse(asm.phase_ids)

        asm.phase_ids.add("c")
        self.assertSequenceEqual(list(asm.phase_ids), ["c"])
        self.assertTrue(asm.phase_ids)
        self.assertFalse("b" in asm.phase_ids)
        self.assertTrue("c" in asm.phase_ids)

    def test_asm_anchors(self):
        asm = bgexvar.Assembly(1000, None, biograph.Sequence("GATTACA"), 5)
        self.assertEqual(asm.seq, "GATTACA")
        self.assertEqual(asm.assembly_id, 5)
        self.assertEqual(asm.left_offset, 1000)
        self.assertEqual(asm.right_offset, None)
        asm.left_offset = None
        asm.right_offset = 1050
        self.assertEqual(asm.left_offset, None)
        self.assertEqual(asm.right_offset, 1050)

    def test_half_anchored_sort(self):
        asm1 = bgexvar.Assembly(999, 1002, biograph.Sequence("GATTACA"), 1)
        asm2 = bgexvar.Assembly(1000, None, biograph.Sequence("GATTACA"), 2)
        asm3 = bgexvar.Assembly(None, 1001, biograph.Sequence("GATTACA"), 3)
        asm4 = bgexvar.Assembly(1002, 1003, biograph.Sequence("GATTACA"), 4)
        for perm in itertools.permutations([asm1, asm2, asm3, asm4]):
            ids = [a.assembly_id for a in sorted(perm)]
            self.assertEqual(ids, [1, 2, 3, 4])

if __name__ == '__main__':
    unittest.main(verbosity=2)
