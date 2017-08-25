#!/usr/bin/env python
import multiprocessing
import os
import sys
from subprocess import Popen, PIPE
from collections import defaultdict

def process_line(tokens):
    path1 = tokens[0].split(':')[0]
    path2 = tokens[1].split(':')[0]
    if path1 == path2: return (tokens, 1, 1)
    exelist = ['benchmark/bin/TMalign', '{}'.format(path1), '{}'.format(path2)]
    proc = Popen(exelist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    actout, acterr = proc.communicate(input = '')
    assert(acterr.strip() == '')
    for outline in actout.split('\n'):
        if outline.startswith('TM-score=') and outline.endswith('(if normalized by length of Chain_1)'):
            tmscore1 = float(outline.split('=')[1].split('(')[0].strip())
        if outline.startswith('TM-score=') and outline.endswith('(if normalized by length of Chain_2)'):
            tmscore2 = float(outline.split('=')[1].split('(')[0].strip())
    return (tokens, tmscore1, tmscore2)

if __name__ == '__main__':
    intokens = sorted([tuple(sorted([line.strip().split()[0], line.strip().split()[1]])) for line in sys.stdin.readlines() if line.strip() != ''])
    if os.path.exists(sys.argv[1]):
        with open(sys.argv[1], 'r') as dbfile:
            dbtokens = set([tuple(sorted([line.strip().split()[0], line.strip().split()[1]])) for line in dbfile.readlines() if line.strip() != ''])
        intokens = [intoken for intoken in intokens if intoken not in dbtokens]
    else:
        with open(sys.argv[1], 'w') as dbfile: pass
    print('Will add {} new entries to the db'.format(len(intokens)))
    results = multiprocessing.Pool().map(process_line, intokens)
    
    os.system('cp "{}" "{}.tmp"'.format(sys.argv[1], sys.argv[1]))
    with open(sys.argv[1]+'.tmp', 'a') as dbfile:    
        for tokens, tmscore1, tmscore2 in results:
            dbfile.write('{}\t{}\t{}\t{}\n'.format(tokens[0], tokens[1], tmscore1, tmscore2))
    os.system('cp "{}.tmp" "{}"'.format(sys.argv[1], sys.argv[1]))

