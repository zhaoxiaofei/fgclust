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
    clusters.add(inner)
    if innfam != outfam:
        n_corrupted_members += 1
        corrupted_clusters.add(inner)
print('Pfam-results: (number-of-clusters, number-of-corrupted-clusters) = ({}, {}) and (number-of-clusters, number-of-corrupted-cluster-members) = ({}, {})'
      .format(len(clusters), len(corrupted_clusters), len(clusters), n_corrupted_members))
#for clu in corrupted_clusters:
#    sys.stdout.write('{}\t'.format(clu))
#sys.stdout.write('\n')

