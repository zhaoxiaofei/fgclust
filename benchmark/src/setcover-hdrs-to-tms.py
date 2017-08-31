#!/usr/bin/env python
import multiprocessing
import sys
from subprocess import Popen, PIPE
from collections import defaultdict

if __name__ == '__main__':
    covering_covered_to_tmscore = {}
    with open(sys.argv[1]) as dbfile:
        lines = [line.strip() for line in dbfile.readlines() if line.strip() != '']
    for line in lines: 
        print(line)
        covering_covered_to_tmscore[(line.split()[0], line.split()[1])] = line.split()[3]
        covering_covered_to_tmscore[(line.split()[1], line.split()[0])] = line.split()[2]
    inner_to_outer_tms_sim_list = defaultdict(list)
    for line in sys.stdin.readlines():
        tokens = line.strip().split()
        tmscore = covering_covered_to_tmscore[(tokens[0], tokens[1])]
        inner_to_outer_tms_sim_list[tokens[0]].append((tokens[1], tmscore, (tokens[2] if len(tokens) > 2 else 'NA')))
    for inner, outer_tms_sim_list in inner_to_outer_tms_sim_list.items():
        outers   = [x[0] for x in outer_tms_sim_list]
        tmscores = [float(x[1]) for x in outer_tms_sim_list]
        sims     = [x[2] for x in outer_tms_sim_list]
        print('names:    min avg {}  {}'.format(inner, outers))
        print('tmscores: {}  {}  N/A {}'.format(min(tmscores), sum(tmscores) / len(tmscores), tmscores))
        print('sims:     N/A N/A N/A {}'.format(sims))

