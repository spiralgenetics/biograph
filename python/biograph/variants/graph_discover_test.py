# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

class ReadCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    @staticmethod
    def only_anchored_vars(asms):
        for a in asms:
            if a.left_offset and a.right_offset and not a.matches_reference:
                yield a

    @staticmethod
    def dedup(asms):
        prev = None
        for a in asms:
            if prev:
                if prev.left_offset == a.left_offset and prev.right_offset == a.right_offset and prev.seq == a.seq:
                    continue
                yield prev
            prev = a
        if prev:
            yield prev

    def trim_and_update(self, asms):
        asms = bgexvar.graph_trim_ref(asms, self.ref, "lambda")
        asms = self.dedup(asms)
        asms = bgexvar.update_rc_seqset_entries(asms, self.seqset, enable_self_test=True)
        return list(asms)

    def test_graph_discover(self):
        asms = bgexvar.make_ref_assemblies(self.ref, "lambda", 0, self.ref.scaffold_lens["lambda"])
        asms = self.trim_and_update(asms)
        asms = bgexvar.discover_branch(asms, self.rm, tag="BRANCH")
        asms = self.trim_and_update(asms)
        asms = bgexvar.discover_push_to_pair(asms, self.rm,
                                             tag="PUSH_TO_PAIR", discover_tags=["BRANCH"])
        asms = self.trim_and_update(asms)
        asms = list(self.only_anchored_vars(asms))

        if len(asms) != 8:
            print("Final list of asms:")
            for a in asms:
                print(a)

        self.assertEqual(len(asms), 8)

if __name__ == '__main__':
    unittest.main(verbosity=2)
