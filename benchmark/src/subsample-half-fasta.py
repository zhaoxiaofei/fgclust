#!/usr/bin/env python
import sys
import random

def fasta_iter(file):
    hdr = seq = None
    for line in file:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            if hdr: yield (hdr, seq)
            hdr = line
            seq = ''
        else: 
            seq += line
    yield (hdr, seq)

random.seed(-1)
for hdr, seq in fasta_iter(sys.stdin):
    if random.randint(0, 1):
        print('{}\n{}'.format(hdr, seq))
    
