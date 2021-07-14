# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import os.path
import tempfile

import biograph
import biograph.variants as bgexvar

def show_progress(new_progress):
    print(f"...{new_progress*100:.2f}")

class ReadCovTestCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bg = biograph.BioGraph("datasets/lambdaToyData/benchmark/father_lambda.bg")
        cls.seqset = cls.bg.seqset
        cls.rm = cls.bg.open_readmap()
        cls.ref = biograph.Reference("datasets/lambdaToyData/benchmark/ref_lambda")

    def test_create_tracer_without_run(self):
        rmap = bgexvar.RefMap.generate_from(self.seqset, self.ref, show_progress)
        disc = bgexvar.ParallelDiscover(self.rm, self.ref, rmap)
        disc.add_entire_reference()
        disc = None

    def test_exception_propagation(self):
        rmap = bgexvar.RefMap.generate_from(self.seqset, self.ref, show_progress)
        disc = bgexvar.ParallelDiscover(self.rm, self.ref, rmap)
        disc.add_entire_reference()
        def output_asm(scaffold_name, asm):
            raise RuntimeError("output_asm")
        with self.assertRaisesRegex(RuntimeError, "output_asm"):
            disc.assemble(output_asm, show_progress)

    def test_discover_parallel(self):
        rmap = bgexvar.RefMap.generate_from(self.seqset, self.ref, show_progress)
        disc = bgexvar.ParallelDiscover(self.rm, self.ref, rmap)
        disc.add_entire_reference()
        asms = []
        disc.assemble(lambda scaffold_name, asm: asms.append((scaffold_name, asm)), show_progress)
        self.assertGreater(len(asms), 5)
        for scaffold_name, asm in asms:
            self.assertEqual(scaffold_name, "lambda")
            self.assertIsInstance(asm, bgexvar.Assembly)

    def test_ref_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rmap_file = tmpdir + "/refmap_file"
            bgexvar.RefMap.generate_and_save(self.seqset, self.ref, rmap_file, show_progress)
            self.assertTrue(os.path.exists(rmap_file))
            rmap = bgexvar.RefMap.load(self.seqset, self.ref, rmap_file)
        disc = bgexvar.ParallelDiscover(self.rm, self.ref, rmap)
        disc.add_entire_reference()
        asms = []
        disc.assemble(lambda scaffold_name, asm: asms.append((scaffold_name, asm)), show_progress)
        self.assertGreater(len(asms), 5)
        for scaffold_name, asm in asms:
            self.assertEqual(scaffold_name, "lambda")
            self.assertIsInstance(asm, bgexvar.Assembly)

if __name__ == '__main__':
    unittest.main(verbosity=2)
