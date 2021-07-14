#!/usr/bin/env python3
from __future__ import print_function
import random
import time
from biograph import BioGraph, Reference


print("Loading biograph...")
mybg = BioGraph('/dev/shm/HG002-30x.bg')

print("Loading reference...")
ref = Reference('/dev/shm/hs37d5')

chromosomes = []
for i in range(22):
    chromosomes.append(str(i+1))

ranges = []
seq_size = 1000

print("Generating reference ranges of size", seq_size)

# Keep trying until we have enough valid ranges
while len(ranges) < 1000:
    chrom = random.choice(chromosomes)
    begin = random.randint(seq_size, ref.scaffold_lens[chrom] - seq_size)
    end = begin + seq_size

    try:
        ranges.append(ref.make_range(chrom, begin, end))
    except RuntimeError:
        print('.', end='')
        continue

print("\nSearching for", len(ranges), "ranges")

now = time.time()
for r in ranges:
    mybg.seq_coverage(r.sequence)
done = time.time()

for r in ranges:
    print(str(r))

print("Elapsed time:", done - now)

# equivalent test for samtools:
# save the above ranges to ranges.txt, and change the format to chr:begin-end
# time for r in `cat ranges.txt`; do samtools depth -r $r /dev/shm/HG002.hs37d5.30x.cram --reference /dev/shm/hs37d5/source.fasta > /dev/null; done
