#!/usr/bin/env python
import sys
for line in sys.stdin:
    tokens = line.strip().split()
    if tokens[0] == 'H':
        print('{}\t{}\t{}'.format(tokens[9], tokens[8], tokens[3]))
    elif tokens[0] == 'C':
        print('{}\t{}\t{}'.format(tokens[8], tokens[8], 100))

