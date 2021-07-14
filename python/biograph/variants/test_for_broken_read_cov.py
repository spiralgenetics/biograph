# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import biograph
import biograph.variants as bgexvar

def vcf_assembly(pos, ref, alt, asm_id):
    pos = int(pos)-1
    if ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return bgexvar.Assembly(pos, pos + len(ref), alt, asm_id)

# Raises an error that F0705 21:53:28.427729 102440 read_cov.cpp:28] Check failed: !m_interior
asms = [vcf_assembly('17125234', 'A', 'G', 17125234),
        vcf_assembly('17125346', 'T', 'C', 17125346),
        vcf_assembly('17125438', 'A', 'C', 17125438),
        vcf_assembly('17125626', 'A', 'C', 17125626)]
bg = biograph.BioGraph("/share/datasets/HG002/HG002-NA24385-50x.bg/", biograph.CacheStrategy.MMAP)
seqset = bg.seqset
rm = bg.open_readmap()
ref = biograph.Reference("/reference/hs37d5")

for pos,asm in enumerate(asms):
    print(pos)
    pc = list(bgexvar.generate_read_cov(rm, ref, "1", [asm]))

