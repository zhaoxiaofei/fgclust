#!/usr/bin/env python
import sys
for line in sys.stdin:
    if line.startswith('#=GF ID'):
        famname = line.strip().split()[2]
    elif line.startswith('#') or line.startswith('//') or line.strip() == '': pass
    else:
        assert(len(line.strip().split()) == 2), '{} is not valid'.format(line)
        faahdr, faaseq = tuple(line.strip().split())
        print('>{}@{}\n{}'.format(famname, faahdr, faaseq.replace('.', '').replace('-', '')))

