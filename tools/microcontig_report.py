#!/usr/bin/env python3
'''
Extract timings from a pcmp console log.
'''

import fileinput
import re

print('microcontig', 'chrom', 'start', 'end', 'seconds', 'variants', 'variants/s', 'bases/s', 'bases', sep='\t')

for line in fileinput.input():

    # 2020-05-12 10:07:21,782 [INFO] 1:329159-331176 done in 5.47s (2 v, 0 v/s, 369 b/s)
    match = re.search(r'^.* \[INFO\] ((.+):(\d+)-(\d+)) done in (.*)s \((\d+) v, (\d+) v/s, (\d+) b/s\)$', line)
    if not match:
        continue

    print('\t'.join(match.groups()), int(match.group(4)) - int(match.group(3)), sep='\t')
