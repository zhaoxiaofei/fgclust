#!/usr/bin/env python
import multiprocessing
import sys
from subprocess import Popen, PIPE
from collections import defaultdict

def process_line(line):
    tokens = line.split()
    path1 = tokens[0].split(':')[0]
    path2 = tokens[1].split(':')[0]
    if path1 == path2: return (tokens, 1)
    exelist = ['benchmark/bin/TMalign', '{}'.format(path1), '{}'.format(path2)]
    proc = Popen(exelist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    actout, acterr = proc.communicate(input = '')
    assert(acterr.strip() == '')
    for outline in actout.split('\n'):
        if outline.startswith('TM-score=') and outline.endswith('(if normalized by length of Chain_2)'):
            tmscore = float(outline.split('=')[1].split('(')[0].strip())
    return (tokens, tmscore)

if __name__ == '__main__':
    lines = [line.strip() for line in sys.stdin.readlines() if line.strip() != '']
    results = multiprocessing.Pool().map(process_line, lines)
    inner_to_outer_tms_sim_list = defaultdict(list)
    for tokens, tmscore in results:
        inner_to_outer_tms_sim_list[tokens[0]].append((tokens[1], tmscore, (tokens[2] if len(tokens) > 2 else 'NA')))
    for inner, outer_tms_sim_list in inner_to_outer_tms_sim_list.items():
        outers   = [x[0] for x in outer_tms_sim_list]
        tmscores = [x[1] for x in outer_tms_sim_list]
        sims     = [x[2] for x in outer_tms_sim_list]
        print('names:    min avg {}  {}'.format(inner, outers))
        print('tmscores: {}  {}  N/A {}'.format(min(tmscores), sum(tmscores) / len(tmscores), tmscores))
        print('sims:     N/A N/A N/A {}'.format(sims))

