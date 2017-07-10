#!/usr/bin/env python
import sys
n_corrupted_members = 0
corrupted_clusters = set([])
clusters = set([])
for line in sys.stdin:
    toks = line.strip().split()
    inner, outer = (toks[0], toks[1])
    innfam = inner.split('@')[0]
    outfam = outer.split('@')[0]
    clusters.add(innfam)
    if innfam != outfam:
        n_corrupted_members += 1
        corrupted_clusters.add(innfam)
print('Pfam-results:{} clusters are corrupted and {} members are corrupted out of {} clusters'.format(len(corrupted_clusters), n_corrupted_members, len(clusters)))

