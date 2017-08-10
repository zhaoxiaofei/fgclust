#!/usr/bin/env python
from __future__ import division
import math
import sys

def process_tmsfname(tmsfname):
    ret = [0] * 10
    with open(tmsfname) as tmsfile:
        for line in tmsfile:
            line = line.strip()
            if line.startswith('tmscore'):
                tms10 = math.ceil(float(line.split()[1]) * 10) - 1
                assert 0 <= tms10 and tms10 <= 9, 'tms10={}'.format(tms10)
                ret[int(tms10)] += 1
    for i in range(1, 10, 1):
        ret[i] = ret[i-1] + ret[i]
    return ret

mine     = process_tmsfname(sys.argv[1])
if len(sys.argv) > 2:
    linclust = process_tmsfname(sys.argv[2])
    cdhit    = process_tmsfname(sys.argv[3])
    quaclust = process_tmsfname(sys.argv[4])
    kclust   = process_tmsfname(sys.argv[5])

#print('TM score & mine & linclust &cd-hit')
#for i in range(10):
#    print('{} & {} & {} & {}'.format(round((i+1)/10, 1), mine[i], linclust[i], cdhit[i]))

sim = 0
if '-50.' in sys.argv[1]:
    sim = 50
if '-70.' in sys.argv[1]:
    sim = 70
if '-90.' in sys.argv[1]:
    sim = 90
assert sim > 0

print('TM score  & {} \\\\'     .format(     ' & '.join([str((i+1)/10.0)  for i in range(10)])))
print('FgClust   & {} & {} \\\\'.format(sim, ' & '.join([str(mine[i])     for i in range(10)])))
if len(sys.argv) > 2:
    print('Linclust  & {} & {} \\\\'.format(sim, ' & '.join([str(linclust[i]) for i in range(10)])))
    print('CD-HIT    & {} & {} \\\\'.format(sim, ' & '.join([str(cdhit[i])    for i in range(10)])))
    print('MMseqs2   & {} & {} \\\\'.format(sim, ' & '.join([str(quaclust[i]) for i in range(10)])))
    print('kClust    & {} & {} \\\\'.format(sim, ' & '.join([str(kclust[i])   for i in range(10)])))

