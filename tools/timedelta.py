#!/usr/bin/env python3.6
'''
Compute the delta between contig processing start/end times on the pcmp console log.

Outputs h:m:s contig base_count
'''

import fileinput
from datetime import datetime

times = dict()

# 2019-07-18 20:25:36,905
fmt = '%Y-%m-%d %H:%M:%S'
for line in fileinput.input():
    try:
        d, hmsm, _, sf, contig = line.split()
    except ValueError:
        continue
    hms, ms = hmsm.split(',')

    thetime = f'{d} {hms}'

    if contig in times:
        start = times[contig]
        try:
            chrom, rng = contig.rsplit(':', 1)
        except ValueError:
            continue
        b, e = rng.split('-')
        print(datetime.strptime(thetime, fmt) - datetime.strptime(start, fmt), contig, int(e) - int(b))
    else:
        times[contig] = thetime
