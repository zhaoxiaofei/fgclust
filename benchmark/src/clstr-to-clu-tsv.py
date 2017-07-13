#!/usr/bin/env python
import sys

class Clu:
    def __init__(self):
        self.inner = None
        self.outer_iden_s = []
    def printout(self):
        print('{}\t{}\t{}'.format(self.inner, self.inner, 100))
        for outer, iden in sorted(self.outer_iden_s):
            print('{}\t{}\t{}'.format(self.inner, outer, iden))

clu = None
for line in sys.stdin:
    line = line.strip()
    if not line: continue
    if line.startswith('>Cluster'):
        if clu: clu.printout()
        clu = Clu()
    elif line.endswith('*'):
        clu.inner = line.split()[2].lstrip('>').rstrip('.')
    elif line.endswith('%'):
        clu.outer_iden_s.append((line.split()[2].lstrip('>').rstrip('.'), float(line.split()[4].rstrip('%').lstrip('+/').lstrip('-/'))))
    else: raise RuntimeError()
if clu: clu.printout()


