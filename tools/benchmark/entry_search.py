#!/usr/bin/env python3
from __future__ import print_function
import random
import time
from biograph import BioGraph

def random_kmer(size=20):
    ret = []
    for i in range(size):
        ret.append(random.choice(['A','C','G','T']))
    return ''.join(ret)

print("Generating kmers...")
kmers = []
for i in range(1000000):
    kmers.append(random_kmer())

print("Loading biograph...")
mybg = BioGraph('HG002-30x.bg')

print("Searching for", len(kmers), "kmers")

now = time.time()
for kmer in kmers:
    mybg.entry_search(kmer)
done = time.time()

print("Elapsed time:", done - now)
