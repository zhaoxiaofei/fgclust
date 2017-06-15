#!/usr/bin/env python
from __future__ import division
import sys

def process_tmsfname(tmsfname, arr):
    with open(tmsfname) as tmsfile:
        for line in tmsfile:
            line = line.strip()
            if line.startswith('tmscore'):
                tms10 = int(round(float(line.split()[1]) * 10, 0))
                arr[tms10] += 1

mine = [0] * 11
linclust = [0] * 11
cdhit = [0] * 11
process_tmsfname(sys.argv[1], mine)
process_tmsfname(sys.argv[2], linclust)
process_tmsfname(sys.argv[3], cdhit)

for i in range(11):
    print('{},{},{},{}'.format(round(i/10, 1), mine[i], linclust[i], cdhit[i]))

