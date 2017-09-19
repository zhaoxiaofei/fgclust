#!/usr/bin/env python
from __future__ import division
import math
import sys

def process_tmsfname(tmsfname):
    ret = [0] * 10
    nclusters = 0
    with open(tmsfname) as tmsfile:
        for line in tmsfile:
            line = line.strip()
            if line.startswith('tmscore'):
                nclusters += 1
                for tok in line.split()[4:]:
                    tms10 = math.ceil(float(tok.lstrip('[').rstrip(']').rstrip(',')) * 10) - 1
                    assert 0 <= tms10 and tms10 <= 9, 'tms10={}'.format(tms10)
                    ret[int(tms10)] += 1
    for i in range(1, 10, 1):
        ret[i] = ret[i-1] + ret[i]
    return (nclusters, ret)

nclu, mine     = process_tmsfname(sys.argv[1])

sim = 0
if '-50.' in sys.argv[1]:
    sim = 50
if '-70.' in sys.argv[1]:
    sim = 70
if '-90.' in sys.argv[1]:
    sim = 90
assert sim > 0

progname = None
if 'hdrsetcover' in sys.argv[1]:
    progname = 'FgClust  '
if 'linclust' in sys.argv[1]:
    progname = 'Linclust '
if 'quaclust' in sys.argv[1]:
    progname = 'MMSeqs2  '
if 'cdhit' in sys.argv[1]:
    progname = 'CD-HIT   '
if 'kclust' in sys.argv[1]:
    progname = 'kClust   '
assert progname

print('TM-score & sim & nclusters & {} \\\\'.format(' & '.join([str((i+1)/10.0)  for i in range(1, 9)])))
print('{} & {} & {} & {} \\\\'.format(progname, sim, nclu, ' & '.join([str(mine[i])     for i in range(1, 9)])))

