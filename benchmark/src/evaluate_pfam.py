import sys

from collections import defaultdict

fam = ''
memtofam = {}
for line in sys.stdin:
    line = line.strip()
    if not line: continue
    if line[0] != '#':
        memtofam[line.split()[0]] = fam
    elif line.startswith('#=GF ID'):
        fam = line.split()[2].strip()
    elif line.startswith('#=GF= AC'):
        fam = fam + '|' + line.split()[2].strip()

corrupted = 0
nseqs = 0
visited = set([])
with open(sys.argv[1]) as faa:
    for line in faa:
        line = line.strip()
        if not line: continue
        if line[0] != '>': continue
        ids = line[1:].split()
        iscorrupted = 0
        for id in ids:
            assert(id not in visited), 'fastaID {} is in two clusters'.format(id)
            visited.add(id)
            if memtofam[id] != memtofam[ids[0]]:
                iscorrupted += 1
        nseqs += len(ids)
        if iscorrupted:
            print('!!!The cluster {} with size {} has {} corrupted members'.format(ids[0], len(ids), iscorrupted))
            corrupted += 1
        else:
            pass
            #print('   The cluster {} with size {} is all in {}'.format(ids[0], len(ids), memtofam[ids[0]]))
print('{} Clusters are corrupted after clustering {} seqs!'.format(corrupted, nseqs))
