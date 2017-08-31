import multiprocessing
import sys
from collections import defaultdict
from dateutil import parser

aatable = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def process_entfname(entfname):
    entfname = entfname.strip()
    faaid_to_seq = defaultdict(list)
    if len(entfname) < 1: return defaultdict(list)
    mod_to_ori = defaultdict()
    with open(entfname) as entfile:
        for line in entfile:
            tokens = line.strip().split()
            try:
                if tokens[0] == 'HEADER' and parser.parse(line[50:59]) > parser.parse('01-JAN-17'): return defaultdict(list)
            except:
                sys.stderr.write('The following line has invalid date: {}\n'.format(line))
            if tokens[0] != 'MODRES': continue
            mod_to_ori[tokens[2]] = tokens[5]
    with open(entfname) as entfile:
        for line in entfile:
            tokens = line.strip().split()
            if len(tokens) < 5 or tokens[0] != 'SEQRES': continue
            faaid = entfname + ':' + tokens[2]
            for aa in tokens[4:]:
                if aa in aatable: faaid_to_seq[faaid].append(aatable[aa])
                elif aa in mod_to_ori and mod_to_ori[aa] in aatable: faaid_to_seq[faaid].append(aatable[mod_to_ori[aa]])
                else: faaid_to_seq[faaid].append('X')
    return faaid_to_seq

if __name__ == '__main__':
    entfnames = [entfname for entfname in sys.stdin]
    #faaid_to_seq_s = map(process_entfname, entfnames)
    faaid_to_seq_s = multiprocessing.Pool().map(process_entfname, entfnames)
    for faaid_to_seq in faaid_to_seq_s:
        if len(faaid_to_seq) != 1: continue
        for faaid, seq in faaid_to_seq.items():
            seq = ''.join(seq)
            if len(seq.replace('X', '')) >= len(seq) * 0.9 and len(seq.replace('X', '')) >= 11: print('>{}\n{}'.format(faaid, seq))

