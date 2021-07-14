# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph

class BioGraphTestCases(unittest.TestCase):

    def test_open_without_cache_strategy(self):
        bg = biograph.BioGraph("golden/e_coli_10000snp.bg")
        self.check_single_metadata(bg)

    def test_open_with_mmap(self):
        bg = biograph.BioGraph(
            "golden/e_coli_10000snp.bg", biograph.CacheStrategy.MMAP)
        self.check_single_metadata(bg)

    def test_open_with_mmapcache(self):
        bg = biograph.BioGraph(
            "golden/e_coli_merged.bg", biograph.CacheStrategy.MMAPCACHE)
        self.check_merged_metadata(bg)

    def test_open_with_ram(self):
        bg = biograph.BioGraph(
            "golden/e_coli_merged.bg", biograph.CacheStrategy.RAM)
        self.check_merged_metadata(bg)

    def check_single_metadata(self, bg):
        md = bg.metadata
        self.assertEqual(md.accession_id, "test_accession_id")
        self.assertEqual(
            md.samples, {"test_accession_id": "87f0fd61143e2441d95c8b8b3ad279473f84be4b"})
        self.assertEqual(md.version, "3.1.1")
        self.assertEqual(
            md.biograph_id, "133cb8e6-f355-9440-38d0-8f9aed4a9ca7")
        rm = bg.open_readmap("")
        self.assertEqual(repr(rm), repr(bg.open_readmap("test_accession_id")))
        self.assertEqual(
            repr(rm), repr(bg.open_readmap("87f0fd61143e2441d95c8b8b3ad279473f84be4b")))

    def check_merged_metadata(self, bg):
        md = bg.metadata
        self.assertEqual(md.accession_id, "merged_id")
        self.assertEqual(md.samples, {
            "test_accession_id": "807ee81ae2650609949f68f4786c9a7dc480f406",
            "e_coli_test": "f21180d242113e56c5bfd645d2614cf16fdde4ed",
        })
        self.assertEqual(md.version, "3.1.1")
        self.assertEqual(
            md.biograph_id, "cd14050f-fcd7-5538-5ab7-f1b8b4051c1c")
        with self.assertRaises(RuntimeError):
            bg.open_readmap("")

        rms = dict()
        for accession_id, uuid in md.samples.items():
            rms[uuid] = bg.open_readmap(accession_id)

        for uuid in md.samples.values():
            self.assertEqual(repr(bg.open_readmap(uuid)), repr(rms[uuid]))

if __name__ == '__main__':
    unittest.main(verbosity=2)
